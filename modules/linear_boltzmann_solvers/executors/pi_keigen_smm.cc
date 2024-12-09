// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/executors/pi_keigen_smm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_pwlc_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_vecops.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/math/parallel_vector/ghosted_parallel_stl_vector.h"
#include "framework/object_factory.h"
#include "framework/utils/timer.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <memory>
#include <numeric>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, PowerIterationKEigenSMM);

InputParameters
PowerIterationKEigenSMM::GetInputParameters()
{
  InputParameters params = PowerIterationKEigen::GetInputParameters();

  params.SetGeneralDescription(
    "Implementation of a second-moment method based k-eigenvalue solver.");
  params.SetDocGroup("LBSExecutors");
  params.ChangeExistingParamToOptional("name", "PowerIterationKEigenSMM");

  // Diffusion k-eigenvalue options
  params.AddOptionalParameter("accel_pi_max_its",
                              50,
                              "The maximum number of power iterations allowed for the "
                              "second-moment method diffusion solver.");
  params.AddOptionalParameter("accel_pi_k_tol",
                              1.0e-10,
                              "The k-eigenvalue tolerance for the second-moment method "
                              "diffusion solver.");
  params.AddOptionalParameter(
    "accel_pi_verbose", false, "A flag for second-moment method diffusion solver verbosity.");

  params.ConstrainParameterRange("accel_pi_k_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("accel_pi_max_its", AllowableRangeLowLimit::New(0));

  // Diffusion linear solve options
  params.AddOptionalParameter("diff_sdm",
                              "pwlc",
                              "The spatial discretization method for the second-moment method "
                              "diffusion system.");
  params.AddOptionalParameter(
    "diff_l_abs_tol", 1.0e-10, "The absolute residual tolerance of the diffusion accelerator.");
  params.AddOptionalParameter(
    "diff_l_max_its", 100, "The maximum allowable iterations for the diffusion accelerator.");
  params.AddOptionalParameter(
    "diff_petsc_options", "", "PETSc options to pass to the diffusion accelerator.");
  params.AddOptionalParameter("diff_verbose", false, "A flag for diffusion accelerator verbosity.");

  params.ConstrainParameterRange("diff_sdm", AllowableRangeList::New({"pwlc", "pwld"}));
  params.ConstrainParameterRange("diff_l_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("diff_l_max_its", AllowableRangeLowLimit::New(0));
  return params;
}

std::shared_ptr<PowerIterationKEigenSMM>
PowerIterationKEigenSMM::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<PowerIterationKEigenSMM>("lbs::PowerIterationKEigenSMM", params);
}

PowerIterationKEigenSMM::PowerIterationKEigenSMM(const InputParameters& params)
  : PowerIterationKEigen(params),
    dimension_(0),
    psi_new_local_(lbs_solver_->PsiNewLocal()),
    accel_pi_max_its_(params.GetParamValue<unsigned int>("accel_pi_max_its")),
    accel_pi_k_tol_(params.GetParamValue<double>("accel_pi_k_tol")),
    accel_pi_verbose_(params.GetParamValue<bool>("accel_pi_verbose")),
    diffusion_sdm_name_(params.GetParamValue<std::string>("diff_sdm")),
    diffusion_l_max_its_(params.GetParamValue<int>("diff_l_max_its")),
    diffusion_l_abs_tol_(params.GetParamValue<double>("diff_l_abs_tol")),
    diffusion_petsc_options_(params.GetParamValue<std::string>("diff_petsc_options")),
    diffusion_verbose_(params.GetParamValue<bool>("diff_verbose"))
{
  ghosts_required_ = diffusion_sdm_name_ == "pwlc";
  lbs_solver_->Options().save_angular_flux = true;
  if (lbs_solver_->Groupsets().size() != 1)
    throw std::logic_error("The SMM k-eigenvalue executor is only implemented for "
                           "problems with a single groupset.");

  // If using the AAH solver with one sweep, a few iterations need to be done
  // to get rid of the junk in the unconverged lagged angular fluxes.  Five
  // sweeps is a guess at how many initial sweeps are necessary.
  auto lbs_solver = std::dynamic_pointer_cast<DiscreteOrdinatesSolver>(lbs_solver_);
  if (lbs_solver->SweepType() == "AAH" and front_gs_.max_iterations == 1)
    throw std::logic_error("The AAH solver is not stable for single-sweep methods due to "
                           "the presence of lagged angular fluxes.  Multiple sweeps are "
                           "allowed, however, the number of sweeps required to get sensible "
                           "results is not well studied and problem dependent.");
}

void
PowerIterationKEigenSMM::Initialize()
{
  PowerIterationKEigen::Initialize();

  // Shorthand information
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& phi_uk_man = lbs_solver_->UnknownManager();
  const auto num_groups = lbs_solver_->NumGroups();
  const auto num_gs_groups = front_gs_.groups.size();

  // Specialized SMM data
  dimension_ = lbs_solver_->Grid().Dimension();
  ComputeAuxiliaryUnitCellMatrices();
  ComputeBoundaryFactors();

  // Create PWLC structure, if needed
  if (diffusion_sdm_name_ == "pwlc")
  {
    pwlc_ptr_ = PieceWiseLinearContinuous::New(grid);
    ghost_info_ = MakePWLDGhostInfo(pwld, phi_uk_man);
  }

  // Create tensor structure
  for (int g = 0; g < num_groups; ++g)
    tensor_uk_man_.AddUnknown(UnknownType::VECTOR_N, dimension_ * dimension_);

  const auto local_size = pwld.GetNumLocalDOFs(tensor_uk_man_);
  const auto global_size = pwld.GetNumGlobalDOFs(tensor_uk_man_);
  const auto ghost_ids = MakePWLDGhostIndices(pwld, tensor_uk_man_);
  tensors_ = std::make_unique<GhostedParallelSTLVector>(
    local_size, global_size, ghost_ids, opensn::mpi_comm);

  // Create diffusion solver
  UnknownManager diff_uk_man;
  diff_uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

  for (const auto& [bid, bc] : lbs_solver_->BoundaryPreferences())
    if ((bc.type == LBSBoundaryType::ISOTROPIC) or (bc.type == LBSBoundaryType::ARBITRARY))
      throw std::logic_error("Only vacuum and reflective boundaries are valid for "
                             "k-eigenvalue problems.");

  const auto bcs = TranslateBCs(lbs_solver_->SweepBoundaries(), false);

  // Create the diffusion materials
  const auto xs_map = PackGroupsetXS(
    lbs_solver_->GetMatID2XSMap(), front_gs_.groups.front().id, front_gs_.groups.back().id);

  // Create the appropriate solver
  log.Log() << "Creating diffusion solver";
  auto& diff_solver = diffusion_solver_;
  if (diffusion_sdm_name_ == "pwld")
  {
    diff_solver = std::make_shared<DiffusionMIPSolver>(std::string(Name() + "_SMM"),
                                                       pwld,
                                                       diff_uk_man,
                                                       bcs,
                                                       xs_map,
                                                       lbs_solver_->GetUnitCellMatrices(),
                                                       true,
                                                       diffusion_verbose_);
  }
  else
  {
    diff_solver = std::make_shared<DiffusionPWLCSolver>(std::string(Name() + "_SMM"),
                                                        *pwlc_ptr_,
                                                        diff_uk_man,
                                                        bcs,
                                                        xs_map,
                                                        lbs_solver_->GetUnitCellMatrices(),
                                                        true,
                                                        diffusion_verbose_);
  }

  diff_solver->options.residual_tolerance = diffusion_l_abs_tol_;
  diff_solver->options.max_iters = diffusion_l_max_its_;
  diff_solver->options.verbose = diffusion_verbose_;
  diff_solver->options.additional_options_string = diffusion_petsc_options_;
  log.Log() << "Done creating diffusion solver.";

  // Initialize the solver
  log.Log() << "Initializing diffusion solver.";
  diff_solver->Initialize();
  mpi_comm.barrier();
  log.Log() << "Done initializing diffusion solver.";

  // Assemble the system
  log.Log() << "Assembling diffusion system.";
  std::vector<double> tmp;
  if (diffusion_sdm_name_ == "pwld")
    tmp.assign(pwld.GetNumLocalDOFs(diff_uk_man), 0.0);
  else
    tmp.assign(pwlc_ptr_->GetNumLocalAndGhostDOFs(diff_uk_man), 0.0);
  AssembleDiffusionBCs();
  diff_solver->AssembleAand_b(tmp);
  log.Log() << "Done assembling diffusion system.";
}

void
PowerIterationKEigenSMM::Execute()
{
  // Start transport power iterations
  double k_eff_ell = k_eff_;
  double k_eff_change = 1.0;

  auto phi_ell = phi_old_local_;

  int nit = 0;
  while (nit < max_iters_)
  {
    // Compute the transport 1/k fission source
    SetLBSFissionSource(phi_ell, false);
    Scale(q_moments_local_, 1.0 / k_eff_);

    // Solve some transport inners
    ags_solver_->Solve();

    std::vector<double> phi0;
    TransferTransportToDiffusion(phi_new_local_, phi0);
    auto phi0_old = phi0;
    auto phi0_m = phi0;

    // Update second-moment method data
    ComputeClosures(psi_new_local_);
    const auto correction = ComputeSourceCorrection();

    // Start diffusion power iterations
    double lambda = k_eff_;
    double lambda_old = k_eff_;
    double F0_old = lbs_solver_->ComputeFissionProduction(phi_new_local_);

    for (int m = 0; m < accel_pi_max_its_; ++m)
    {
      // Define the nodal diffusion 1/lambda fission source.
      // This returns a local PWLD vector which must still be integrated
      // against the diffusion test functions.
      auto Sf = SetNodalDiffusionFissionSource(phi0_m);
      Scale(Sf, 1.0 / lambda);

      // Solve the diffusion inners. This seems quite wasteful, so
      // only a single iteration is taken.
      for (int i = 0; i < 1; ++i)
      {
        // Define the nodal diffusion scattering source.
        // This returns a local PWLD vector which must still be integrated
        // against the diffusion test functions.
        auto Ss = SetNodalDiffusionScatterSource(phi0_old);

        // Assemble the normal terms of the RHS of the diffusion solver,
        // add in the second moment correction, and solve the system
        auto b = AssembleDiffusionRHS(Sf + Ss);
        diffusion_solver_->AddToRHS(b + correction);
        diffusion_solver_->Solve(phi0, true);

        // Bump the new solution to old
        phi0_old = phi0;
      }

      // Compute a new diffusion eigenvalue. This requires mapping the
      // diffusion solution back to the transport discretization.
      std::vector<double> phi;
      TransferDiffusionToTransport(phi0, phi);
      double F0 = lbs_solver_->ComputeFissionProduction(phi);
      lambda = F0 / F0_old * lambda_old;

      // Check for convergence
      double lambda_change = std::fabs(lambda / lambda_old - 1.0);
      if (accel_pi_verbose_)
        log.Log() << "SMM PI iteration " << m << "  lambda " << lambda << "  change "
                  << lambda_change << (lambda_change < accel_pi_k_tol_ ? "  CONVERGED" : "");

      // Bump variables for next iteration
      lambda_old = lambda;
      F0_old = F0;
      phi0_m = phi0;

      if (lambda_change < accel_pi_k_tol_)
        break;

    } // acceleration

    // Project back to transport and check convergence
    k_eff_ = lambda;
    k_eff_change = std::fabs(k_eff_ / k_eff_ell - 1.0);
    k_eff_ell = k_eff_;

    TransferDiffusionToTransport(phi0, phi_new_local_);

    const double production = lbs_solver_->ComputeFissionProduction(phi_new_local_);
    LBSVecOps::ScalePhiVector(*lbs_solver_, PhiSTLOption::PHI_NEW, lambda / production);

    const auto phi_change = CheckScalarFluxConvergence(phi_new_local_, phi_ell);
    LBSVecOps::GSScopedCopyPrimarySTLvectors(
      *lbs_solver_, front_gs_, phi_new_local_, phi_old_local_);
    LBSVecOps::GSScopedCopyPrimarySTLvectors(*lbs_solver_, front_gs_, phi_new_local_, phi_ell);
    ++nit;

    if (lbs_solver_->Options().verbose_outer_iterations)
    {
      std::stringstream ss;
      ss << program_timer.GetTimeString() << " "
         << "  Iteration " << std::setw(5) << nit << "  k_eff " << std::setw(14)
         << std::setprecision(10) << k_eff_ << "  k_eff change " << std::setw(12) << k_eff_change
         << "  phi change " << std::setw(12) << phi_change
         << (k_eff_change < k_tolerance_ ? " CONVERGED" : "");
      log.Log() << ss.str();
    }

    if (k_eff_change < k_tolerance_)
      break;
  }

  // Print summary
  log.Log() << "\n";
  log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff_;
  log.Log() << "        Final change          :        " << std::setprecision(6) << k_eff_change
            << " (Number of Sweeps:" << front_wgs_context_->counter_applications_of_inv_op << ")"
            << "\n";

  if (lbs_solver_->Options().use_precursors)
  {
    lbs_solver_->ComputePrecursors();
    Scale(lbs_solver_->PrecursorsNewLocal(), 1.0 / k_eff_);
  }

  lbs_solver_->UpdateFieldFunctions();
}

void
PowerIterationKEigenSMM::ComputeClosures(const std::vector<std::vector<double>>& psi)
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& transport_views = lbs_solver_->GetCellTransportViews();

  // Create a local tensor vector, set it to zero
  auto local_tensors = tensors_->MakeLocalVector();
  local_tensors.assign(local_tensors.size(), 0.0);

  // Loop over groupsets
  int gs = 0;
  for (const auto& groupset : lbs_solver_->Groupsets())
  {
    const auto& psi_uk_man = groupset.psi_uk_man_;
    const auto& quad = groupset.quadrature;
    const auto num_gs_dirs = quad->omegas.size();

    const auto first_grp = groupset.groups.front().id;
    const auto num_gs_groups = groupset.groups.size();

    // Loop over cells
    for (const auto& cell : grid.local_cells)
    {
      const auto& transport_view = transport_views[cell.local_id];
      const auto& cell_mapping = pwld.GetCellMapping(cell);

      // Compute node-wise, groupset wise tensors
      for (int i = 0; i < transport_view.NumNodes(); ++i)
      {
        for (int gsg = 0; gsg < num_gs_groups; ++gsg)
        {
          const auto g = first_grp + gsg;

          // Get the tensor, set it to zero
          const auto dof_map = pwld.MapDOFLocal(cell, i, tensor_uk_man_, g, 0);
          double* T = &local_tensors[dof_map];

          // Perform the quad integration
          for (int d = 0; d < num_gs_dirs; ++d)
          {
            const auto& omega = quad->omegas[d];
            const auto& wt = quad->weights[d];

            const auto psi_dof = pwld.MapDOFLocal(cell, i, psi_uk_man, d, gsg);
            const auto coeff = wt * psi[gs][psi_dof];

            for (int k = 0; k < dimension_; ++k)
            {
              const auto dim_idx_k = dimension_ > 1 ? k : 2;
              T[k * dimension_ + k] -= coeff / 3.0;

              for (int l = 0; l < dimension_; ++l)
              {
                const auto dim_idx_l = dimension_ > 1 ? l : 2;
                T[k * dimension_ + l] += coeff * omega[dim_idx_k] * omega[dim_idx_l];
              }
            }
          } // for direction d
        }   // for groupset group gsg
      }     // for node i

      // Loop over cell faces
      int f = 0;
      for (const auto& face : cell.faces)
      {
        if (not face.has_neighbor)
        {
          const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
          for (int fi = 0; fi < num_face_nodes; ++fi)
          {
            const auto i = cell_mapping.MapFaceNode(f, fi);
            const auto imap = pwld.MapDOFLocal(cell, i);
            const auto bfac = bndry_factors_[imap][gs];

            // Reset boundary closures
            betas_[imap].assign(lbs_solver_->NumGroups(), 0.0);
            auto& beta = betas_[imap];

            // Compute the boundary closure
            for (int gsg = 0; gsg < num_gs_groups; ++gsg)
            {
              const auto g = first_grp + gsg;

              // Perform the quadrature integration
              for (int d = 0; d < num_gs_dirs; ++d)
              {
                const auto& wt = quad->weights[d];
                const auto& omega = quad->omegas[d];
                const auto mu = std::fabs(omega.Dot(face.normal));

                const auto psi_dof = pwld.MapDOFLocal(cell, i, psi_uk_man, d, gsg);
                beta[g] += wt * (mu - bfac) * psi[gs][psi_dof];
              } // for direction n
            }   // for groupset group gsg
          }     // for face node fi
        }
        ++f;
      } // for face
    }   // for cell
    ++gs;
  } // for groupset

  // Set the ghosted tensor vector with the local tensors, then
  // communicate the ghosts
  tensors_->Set(local_tensors);
  tensors_->CommunicateGhostEntries();
}

std::vector<double>
PowerIterationKEigenSMM::ComputeSourceCorrection() const
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& matid_to_xs_map = lbs_solver_->GetMatID2XSMap();
  const auto& unit_cell_matrices = lbs_solver_->GetUnitCellMatrices();

  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  auto tensors = tensors_->MakeGhostedLocalVector();

  // Create the output vector
  ParallelSTLVector output(
    diff_sd.GetNumLocalDOFs(diff_uk_man), diff_sd.GetNumGlobalDOFs(diff_uk_man), mpi_comm);

  // Build the source
  for (const auto& cell : grid.local_cells)
  {
    const auto& rho = lbs_solver_->DensitiesLocal()[cell.local_id];
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto nodes = cell_mapping.GetNodeLocations();
    const auto num_cell_nodes = cell_mapping.GetNumNodes();
    const auto num_cell_faces = cell.faces.size();

    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& K = K_tensor_matrices_[cell.local_id];

    const auto& xs = matid_to_xs_map.at(cell.material_id);
    const auto& sigma_tr = xs->GetSigmaTransport();

    // Volumetric term
    for (int gsg = 0; gsg < num_gs_groups; ++gsg)
    {
      const auto g = first_grp + gsg;
      const auto& sig_tr = rho * sigma_tr[g];

      // Compute the volumetric correction term to the diffusion source
      // for this diffusion DoF. Unlike in many other areas of the codebase,
      // these tensors are (dim x dim) instead of (3 x 3) so no special
      // treatment is needed for the 1D case when only the z-component is
      // present. The term is given by:
      // grad bi dot div(bj Tj) / sigma_tr

      for (int i = 0; i < num_cell_nodes; ++i)
      {
        const auto imap = diff_sd.MapDOF(cell, i, diff_uk_man, 0, gsg);

        double val = 0.0;
        for (int j = 0; j < num_cell_nodes; ++j)
        {
          // Here, because the tensors are transport-based, the transport
          // discretization is used to obtain the tensor DoF.
          const auto jmap = pwld.MapDOFLocal(cell, j, tensor_uk_man_, g, 0);
          const double* T = &tensors[jmap];

          // This adds the rank-2 tensor contraction (double dot product)
          // of the closure tensor and the tensor stiffness matrix for nodes
          // i and j.
          for (int k = 0; k < dimension_; ++k)
            for (int l = 0; l < dimension_; ++l)
              val += T[k * dimension_ + l] * K(i, j, k, l);
        } // for node j

        output.SetValue(imap, -val / sig_tr, VecOpType::ADD_VALUE);

      } // for node i
    }   // for groupset group gsg

    // Surface terms
    for (int f = 0; f < num_cell_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& normal = face.normal;
      const auto& face_G = fe_values.intS_shapeI_gradshapeJ[f];
      const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);

      // Interior face terms
      if (face.has_neighbor)
      {
        const auto& nbr_cell = grid.cells[face.neighbor_id];
        const auto& nbr_cell_mapping = pwld.GetCellMapping(nbr_cell);
        const auto& nbr_nodes = nbr_cell_mapping.GetNodeLocations();

        for (int gsg = 0; gsg < num_gs_groups; ++gsg)
        {
          const auto g = first_grp + gsg;
          const auto& sig_tr = sigma_tr[g];

          // Contribute the closure tensor jump term. This is used in both
          // PWLD and PWLC discretizations since the closure tensors are
          // discontinuous. The loop over cell nodes i is used because the
          // gradient of all basis functions on the current cell are non-zero
          // on the face. This term is given by:
          // 0.5 (Tj^+ bj^+ - Tj^- bj^-) n dot grad bi^- / sig_tr^-

          for (int i = 0; i < num_cell_nodes; ++i)
          {
            const auto imap = diff_sd.MapDOF(cell, i, diff_uk_man, 0, gsg);

            double val = 0.0;
            for (int fj = 0; fj < num_face_nodes; ++fj)
            {
              // Get the m(-) and p(+) cell node indices for this face node.
              // Here the m(-) indicates the current face (the cell in the
              // direction opposite the face normal), and the p(+) indicates
              // the neighbor cell.
              const auto jm = cell_mapping.MapFaceNode(f, fj);
              const auto jp = MapAssociatedFaceNode(nodes[jm], nbr_nodes);

              // Get the SM tensor DoFs
              const auto jmmap = pwld.MapDOFLocal(cell, jm, tensor_uk_man_, g, 0);
              const auto jpmap =
                grid.IsCellLocal(face.neighbor_id)
                  ? pwld.MapDOFLocal(nbr_cell, jp, tensor_uk_man_, g, 0)
                  : tensors_->MapGhostToLocal(pwld.MapDOF(nbr_cell, jp, tensor_uk_man_, g, 0));

              const double* Tm = &tensors[jmmap];
              const double* Tp = &tensors[jpmap];

              for (int k = 0; k < dimension_; ++k)
              {
                const auto dim_idx_k = dimension_ > 1 ? k : 2;
                for (int l = 0; l < dimension_; ++l)
                {
                  const auto dim_idx_l = dimension_ > 1 ? l : 2;
                  const auto dT = Tm[k * dimension_ + l] - Tp[k * dimension_ + l];
                  val += dT * normal[dim_idx_l] * face_G(jm, i)[dim_idx_k];
                }
              }
            } // for face node fj

            output.SetValue(imap, 0.5 * val / sig_tr, VecOpType::ADD_VALUE);

          } // for node i

          // Contribute the average closure tensor divergence term.
          // This is only used in PWLD discretizations because in PWLC there
          // is no jump in the basis functions across a face. The loop over
          // cell nodes j is used because the gradient of all basis functions
          // on the current cell are non-zero on the face.
          // This term is given by:
          // 0.5 Tj^- (bi^+ - bi^-) grad bj^- dot n / sig_tr

          if (diffusion_sdm_name_ == "pwld")
          {
            for (int fi = 0; fi < num_face_nodes; ++fi)
            {
              // Get the current and neighbor cell node indices
              const auto im = cell_mapping.MapFaceNode(f, fi);
              const auto ip = MapAssociatedFaceNode(nodes[im], nbr_nodes);

              // Get the corresponding current and neighbor cell DoFs
              const auto immap = diff_sd.MapDOF(cell, im, diff_uk_man, 0, gsg);
              const auto ipmap = diff_sd.MapDOF(nbr_cell, ip, diff_uk_man, 0, gsg);

              double val = 0.0;
              for (int j = 0; j < num_cell_nodes; ++j)
              {
                const auto jmap = pwld.MapDOFLocal(cell, j, tensor_uk_man_, g, 0);
                const double* Tm = &tensors[jmap];

                for (int k = 0; k < dimension_; ++k)
                {
                  const auto dim_idx_k = dimension_ > 1 ? k : 2;
                  for (int l = 0; l < dimension_; ++l)
                  {
                    const auto dim_idx_l = dimension_ > 1 ? l : 2;
                    val += Tm[k * dimension_ + l] * face_G(im, j)[dim_idx_l] * normal[dim_idx_k];
                  }
                }
              } // for node j

              val *= 0.5 / sig_tr;
              output.SetValue(immap, -val, VecOpType::ADD_VALUE);
              output.SetValue(ipmap, val, VecOpType::ADD_VALUE);

            } // for face node fi
          }
        } // for groupset group gsg
      }   // if face has neighbor

      // Boundary face terms
      if (not face.has_neighbor)
      {
        BoundaryCondition bc;
        if (diffusion_solver_->BCS().count(face.neighbor_id) == 1)
          bc = diffusion_solver_->BCS().at(face.neighbor_id);

        // Contribute the boundary closure term to the diffusion source term.
        // This term is the sum of the integrated incident and outgoing
        // partial currents minus half the quadrature computed scalar flux.
        if (bc.type == BCType::ROBIN)
        {
          // Skip reflective boundaries
          if (bc.values[1] == 1.0)
            continue;

          const auto& face_M = fe_values.intS_shapeI_shapeJ[f];
          for (int gsg = 0; gsg < num_gs_groups; ++gsg)
          {
            const auto g = first_grp + gsg;
            for (int fi = 0; fi < num_face_nodes; ++fi)
            {
              const auto i = cell_mapping.MapFaceNode(f, fi);
              const auto imap = diff_sd.MapDOF(cell, i, diff_uk_man, 0, gsg);

              double val = 0.0;
              for (int fj = 0; fj < num_face_nodes; ++fj)
              {
                const auto j = cell_mapping.MapFaceNode(f, fj);
                const auto jmap = pwld.MapDOFLocal(cell, j);
                val += betas_.at(jmap)[g] * face_M(i, j);
              } // for face node fj

              output.SetValue(imap, -val, VecOpType::ADD_VALUE);

            } // for face node fi
          }   // for groupset group gsg
        }     // if Robin boundary
      }       // if face is boundary
    }         // for face f
  }           // for cell

  output.Assemble();
  return output.MakeLocalVector();
}

void
PowerIterationKEigenSMM::AssembleDiffusionBCs() const
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& unit_cell_matrices = lbs_solver_->GetUnitCellMatrices();

  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();

  const auto num_gs_groups = front_gs_.groups.size();

  // Loop over cells
  std::vector<int64_t> rows;
  std::vector<int64_t> cols;
  std::vector<double> vals;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto& fe_values = unit_cell_matrices[cell.local_id];

    // Loop over faces
    int f = 0;
    for (const auto& face : cell.faces)
    {
      if (not face.has_neighbor)
      {
        BoundaryCondition bc;
        if (diffusion_solver_->BCS().count(face.neighbor_id))
          bc = diffusion_solver_->BCS().at(face.neighbor_id);

        if (bc.type == BCType::ROBIN)
        {
          // skip reflective bcs
          if (std::fabs(bc.values[0]) < 1.0e-12)
            continue;

          const auto& face_M = fe_values.intS_shapeI_shapeJ[f];

          const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
          for (int fi = 0; fi < num_face_nodes; ++fi)
          {
            const auto i = cell_mapping.MapFaceNode(f, fi);
            const auto imap = diff_sd.MapDOF(cell, i, diff_uk_man, 0, 0);

            for (int fj = 0; fj < num_face_nodes; ++fj)
            {
              const auto j = cell_mapping.MapFaceNode(f, fj);
              const auto jmap = diff_sd.MapDOF(cell, j, diff_uk_man, 0, 0);
              const auto bfac = bndry_factors_.at(pwld.MapDOFLocal(cell, j))[0];

              for (int gsg = 0; gsg < num_gs_groups; ++gsg)
              {
                // This will not work if extended to multiple groupsets because in
                // the current paradigm, one diffusion solver per groupset would be
                // required.
                rows.push_back(imap + gsg);
                cols.push_back(jmap + gsg);
                vals.push_back(bfac * face_M(i, j));
              }
            } // for face node fj
          }   // for face node fi
        }     // if Robin boundary
      }       // if boundary face
      ++f;
    } // for face
  }   // for cell

  // Add the contributions to the diffusion solver matrix
  diffusion_solver_->AddToMatrix(rows, cols, vals);
}

std::vector<double>
PowerIterationKEigenSMM::AssembleDiffusionRHS(const std::vector<double>& q0) const
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& uk_man = lbs_solver_->UnknownManager();
  const auto& unit_cell_matrices = lbs_solver_->GetUnitCellMatrices();

  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  // Clear the diffusion RHS
  const auto num_local_diff_dofs = ghosts_required_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man)
                                                    : diff_sd.GetNumLocalDOFs(diff_uk_man);
  std::vector<double> dummy(num_local_diff_dofs, 0.0);
  diffusion_solver_->Assemble_b(dummy);

  // Create the output vector
  ParallelSTLVector rhs(
    diff_sd.GetNumLocalDOFs(diff_uk_man), diff_sd.GetNumGlobalDOFs(diff_uk_man), mpi_comm);

  // Test the nodal source against diffusion test functions
  for (const auto& cell : grid.local_cells)
  {
    const auto& M = unit_cell_matrices[cell.local_id].intV_shapeI_shapeJ;
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto num_cell_nodes = cell_mapping.GetNumNodes();
    for (int gsg = 0; gsg < num_gs_groups; ++gsg)
    {
      const auto g = first_grp + gsg;
      for (int i = 0; i < num_cell_nodes; ++i)
      {
        const auto imap = diff_sd.MapDOF(cell, i, diff_uk_man, 0, gsg);

        double val = 0.0;
        for (int j = 0; j < num_cell_nodes; ++j)
        {
          const auto jmap = pwld.MapDOFLocal(cell, j, uk_man, 0, g);
          val += M(i, j) * q0[jmap];
        }
        rhs.SetValue(imap, val, VecOpType::ADD_VALUE);
      } // for node i
    }   // for group gsg
  }     // for cell

  rhs.Assemble();
  return rhs.MakeLocalVector();
}

std::vector<double>
PowerIterationKEigenSMM::SetNodalDiffusionFissionSource(const std::vector<double>& phi0) const
{
  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto num_local_diff_dofs = ghosts_required_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man)
                                                    : diff_sd.GetNumLocalDOFs(diff_uk_man);

  if (phi0.size() != num_local_diff_dofs)
    throw std::invalid_argument("Vector size mismatch. The flux used to construct the "
                                "diffusion fission source must have the same size as the "
                                "local diffusion solver.");

  std::vector<double> phi;
  TransferDiffusionToTransport(phi0, phi);

  std::vector<double> Sf(phi.size(), 0.0);
  active_set_source_function_(
    front_gs_, Sf, phi, APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  return Sf;
}

std::vector<double>
PowerIterationKEigenSMM::SetNodalDiffusionScatterSource(const std::vector<double>& phi0) const
{
  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();
  const auto num_local_diff_dofs = ghosts_required_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man)
                                                    : diff_sd.GetNumLocalDOFs(diff_uk_man);

  if (phi0.size() != num_local_diff_dofs)
    throw std::invalid_argument("Vector size mismatch. The flux used to construct the "
                                "diffusion fission source must have the same size as the "
                                "local diffusion solver.");

  std::vector<double> phi;
  TransferDiffusionToTransport(phi0, phi);

  std::vector<double> Ss(phi.size(), 0.0);
  active_set_source_function_(front_gs_,
                              Ss,
                              phi,
                              APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
                                SUPPRESS_WG_SCATTER);
  return Ss;
}

void
PowerIterationKEigenSMM::ComputeAuxiliaryUnitCellMatrices()
{
  const auto& discretization = lbs_solver_->SpatialDiscretization();

  // Spatial weight functions
  struct SpatialWeightFunction
  {
    virtual ~SpatialWeightFunction() = default;
    virtual double operator()(const Vector3& p) const { return 1.0; }
  };

  struct CylindricalWeightFunction : public SpatialWeightFunction
  {
    double operator()(const Vector3& p) const override { return p[0]; }
  };

  struct SphericalWeightFunction : public SpatialWeightFunction
  {
    double operator()(const Vector3& p) const override { return p[2] * p[2]; }
  };

  auto swf = std::make_shared<SpatialWeightFunction>();
  const auto geom_type = lbs_solver_->Options().geometry_type;
  if (geom_type == GeometryType::ONED_SPHERICAL)
    swf = std::make_shared<SphericalWeightFunction>();
  else if (geom_type == GeometryType::TWOD_CYLINDRICAL)
    swf = std::make_shared<CylindricalWeightFunction>();

  // Compute integrals
  const auto num_local_cells = lbs_solver_->Grid().local_cells.size();
  K_tensor_matrices_.resize(num_local_cells);
  for (const auto& cell : lbs_solver_->Grid().local_cells)
  {
    const auto& cell_mapping = discretization.GetCellMapping(cell);
    const auto num_cell_nodes = cell_mapping.GetNumNodes();
    const auto qp_data = cell_mapping.MakeVolumetricFiniteElementData();

    NDArray<double, 4> K(
      std::array<size_t, 4>({num_cell_nodes, num_cell_nodes, dimension_, dimension_}), 0.);
    for (int i = 0; i < num_cell_nodes; ++i)
      for (int j = 0; j < num_cell_nodes; ++j)
      {
        for (int k = 0; k < dimension_; ++k)
          for (int l = 0; l < dimension_; ++l)
            for (const auto& qp : qp_data.GetQuadraturePointIndices())
              K(i, j, k, l) += (*swf)(qp_data.QPointXYZ(qp)) *
                               qp_data.ShapeGrad(i, qp)[dimension_ > 1 ? k : 2] *
                               qp_data.ShapeGrad(j, qp)[dimension_ > 1 ? l : 2] * qp_data.JxW(qp);
      }
    K_tensor_matrices_[cell.local_id] = K;
  }

  opensn::mpi_comm.barrier();
  if (accel_pi_verbose_)
    log.Log() << "Second moment method cell matrices computed.";
}

void
PowerIterationKEigenSMM::ComputeBoundaryFactors()
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto num_groupsets = lbs_solver_->Groupsets().size();

  // Loop over groupsets
  int gs = 0;
  for (const auto& groupset : lbs_solver_->Groupsets())
  {
    const auto& quad = groupset.quadrature;
    const auto num_gs_dirs = quad->omegas.size();
    const auto wt_sum = std::accumulate(quad->weights.begin(), quad->weights.end(), 0.0);

    // Loop over cells
    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = pwld.GetCellMapping(cell);

      // Loop over faces
      int f = 0;
      for (const auto& face : cell.faces)
      {
        const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
        for (int fi = 0; fi < num_face_nodes; ++fi)
        {
          const auto i = cell_mapping.MapFaceNode(f, fi);
          const auto imap = pwld.MapDOFLocal(cell, i);

          if (bndry_factors_.count(imap) == 0)
            bndry_factors_[imap].resize(num_groupsets, 0.0);

          // Compute the diffusion boundary factor, or the quadrature integral of
          // |\Omega \cdot \hat{n}| divided by the sum of the quadrature weights.
          double val = 0.0;
          for (int d = 0; d < num_gs_dirs; ++d)
            val += quad->weights[d] * std::fabs(quad->omegas[d].Dot(face.normal));
          bndry_factors_[imap][gs] = val / wt_sum;
        }
        ++f;
      } // for face
    }   // for cell
    ++gs;
  } // for groupset
}

void
PowerIterationKEigenSMM::TransferTransportToDiffusion(const std::vector<double>& input,
                                                      std::vector<double>& output) const
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();

  const auto& phi_uk_man = lbs_solver_->UnknownManager();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  const auto num_local_diff_dofs = ghosts_required_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man)
                                                    : diff_sd.GetNumLocalDOFs(diff_uk_man);

  // Check the input vector
  if (input.size() != pwld.GetNumLocalDOFs(phi_uk_man))
    throw std::logic_error("Vector size mismatch. Expected an input vector of size " +
                           std::to_string(pwld.GetNumLocalDOFs(phi_uk_man)) + ", but got a " +
                           "vector of size " + std::to_string(input.size()) + ".");

  // If diffusion is PWLC, then nodal averages should be transferred to the PWLC vector,
  // otherwise, the result is a copy
  std::vector<double> mapped_phi;
  if (pwlc_ptr_)
    mapped_phi = ComputeNodallyAveragedPWLDVector(input, pwld, diff_sd, phi_uk_man, ghost_info_);
  else
    mapped_phi = input;

  // Go through the cells and transfer the data to the output vector
  output.assign(num_local_diff_dofs, 0.0);
  for (const auto& cell : grid.local_cells)
  {
    // This is the same for PWLC and PWLD, so only this is needed
    const auto& cell_mapping = pwld.GetCellMapping(cell);

    for (int i = 0; i < cell_mapping.GetNumNodes(); ++i)
    {
      const auto input_dof_map = pwld.MapDOFLocal(cell, i, phi_uk_man, 0, first_grp);
      const auto output_dof_map = diff_sd.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      std::copy_n(&mapped_phi[input_dof_map], num_gs_groups, &output[output_dof_map]);
    }
  }
}

void
PowerIterationKEigenSMM::TransferDiffusionToTransport(const std::vector<double>& input,
                                                      std::vector<double>& output) const
{
  const auto& grid = lbs_solver_->Grid();
  const auto& pwld = lbs_solver_->SpatialDiscretization();
  const auto& diff_sd = diffusion_solver_->SpatialDiscretization();

  const auto& uk_man = lbs_solver_->UnknownManager();
  const auto& diff_uk_man = diffusion_solver_->UnknownStructure();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  const auto num_local_diff_dofs = ghosts_required_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man)
                                                    : diff_sd.GetNumLocalDOFs(diff_uk_man);

  // Check the input vector
  if (input.size() != num_local_diff_dofs)
    throw std::logic_error("Vector size mismatch. Expected an input vector of size " +
                           std::to_string(diff_sd.GetNumLocalDOFs(diff_uk_man)) +
                           ", but got a vector of size " + std::to_string(input.size()) + ".");

  // Go through the cells and transfer data
  output.assign(pwld.GetNumLocalDOFs(uk_man), 0.0);
  for (const auto& cell : grid.local_cells)
  {
    // This is the same for PWLC and PWLD, so only this is needed
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto num_cell_nodes = cell_mapping.GetNumNodes();

    for (int i = 0; i < num_cell_nodes; ++i)
    {
      const auto input_dof_map = diff_sd.MapDOFLocal(cell, i, diff_uk_man, 0, 0);
      const auto output_dof_map = pwld.MapDOFLocal(cell, i, uk_man, 0, first_grp);
      std::copy_n(&input[input_dof_map], num_gs_groups, &output[output_dof_map]);
    }
  }
}

double
PowerIterationKEigenSMM::CheckScalarFluxConvergence(const std::vector<double>& phi_new,
                                                    const std::vector<double>& phi_old)
{
  const auto& grid = lbs_solver_->Grid();
  const auto& sdm = lbs_solver_->SpatialDiscretization();
  const auto& uk_man = lbs_solver_->UnknownManager();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  double local_l2norm = 0.0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto num_cell_nodes = cell_mapping.GetNumNodes();

    for (int i = 0; i < num_cell_nodes; ++i)
      for (int gsg = 0; gsg < num_gs_groups; ++gsg)
      {
        const auto g = first_grp + gsg;
        const auto dof = sdm.MapDOFLocal(cell, i, uk_man, 0, g);

        const double dphi = phi_new.at(dof) - phi_old.at(dof);
        local_l2norm += dphi * dphi;
      }
  }

  double global_l2norm = 0.0;
  mpi_comm.all_reduce(local_l2norm, global_l2norm, mpicpp_lite::op::sum<double>());
  return std::sqrt(global_l2norm);
}

std::vector<double>
PowerIterationKEigenSMM::ComputeNodallyAveragedPWLDVector(const std::vector<double>& pwld_vector,
                                                          const SpatialDiscretization& pwld,
                                                          const SpatialDiscretization& pwlc,
                                                          const UnknownManager& uk_man,
                                                          const GhostInfo& ghost_info)
{
  const auto& ghost_comm = ghost_info.vector_ghost_communicator;
  const auto& pwld_global_to_local_map = ghost_info.ghost_global_to_local_map;

  auto ghosted_pwld_vector = ghost_comm->MakeGhostedVector(pwld_vector);
  ghost_comm->CommunicateGhostEntries(ghosted_pwld_vector);

  const auto& grid = pwld.GetGrid();
  const auto num_local_pwlc_dofs = pwlc.GetNumLocalAndGhostDOFs(uk_man);

  // Map a PWLD vector to a PWLC vector by summing all PWLD DoFs into
  // their associated PWLC DoF. Each time a PWLC DoF is encountered, an
  // associated counter is incremented for later averaging
  std::vector<double> pwlc_wtd_avg(num_local_pwlc_dofs, 0.0);
  std::vector<double> pwlc_vol_sum(num_local_pwlc_dofs, 0.0);
  std::map<int64_t, int64_t> pwlc_global_to_local_map;

  // Add local cell data
  std::set<uint64_t> partition_vertex_ids;
  for (const auto& cell : grid.local_cells)
  {
    // Add PWLD DoFs into a PWLC vector. Each time that a PWLD DoF is
    // added into a PWLC DoF, a counter for that DoF is incremented so
    // that the PWLC DoFs can later be averaged.
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto& vol = cell_mapping.GetCellVolume();
    for (int i = 0; i < cell_mapping.GetNumNodes(); ++i)
      for (int u = 0; u < uk_man.unknowns.size(); ++u)
        for (int c = 0; c < uk_man.unknowns[u].num_components; ++c)
        {
          const auto pwld_dof = pwld.MapDOFLocal(cell, i, uk_man, u, c);
          const auto pwlc_dof = pwlc.MapDOFLocal(cell, i, uk_man, u, c);
          const auto pwlc_global_dof = pwlc.MapDOF(cell, i, uk_man, u, c);

          pwlc_global_to_local_map[pwlc_global_dof] = pwlc_dof;
          pwlc_wtd_avg[pwlc_dof] += pwld_vector[pwld_dof] * vol;
          pwlc_vol_sum[pwlc_dof] += vol;
        }

    // Determine vertices that lie on a partition boundary. This is done
    // to later add PWLD DoFs from non-locally owned cells that correspond to
    // a PWLC DoF
    for (const auto& face : cell.faces)
      if (face.has_neighbor)
        if (not grid.IsCellLocal(face.neighbor_id))
          for (const uint64_t vertex_id : face.vertex_ids)
            partition_vertex_ids.insert(vertex_id);
  } // for local cell

  // Add ghost cell data
  const auto ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  const auto& pvids = partition_vertex_ids;
  for (const uint64_t global_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[global_id];
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto& vol = cell_mapping.GetCellVolume();

    for (int i = 0; i < cell_mapping.GetNumNodes(); ++i)
      if (pvids.find(cell.vertex_ids[i]) != pvids.end())
        for (int u = 0; u < uk_man.unknowns.size(); ++u)
          for (int c = 0; c < uk_man.unknowns[u].num_components; ++c)
          {
            const auto pwld_global_dof = pwld.MapDOF(cell, i, uk_man, u, c);
            const auto pwlc_global_dof = pwlc.MapDOF(cell, i, uk_man, u, c);
            if (pwlc_global_to_local_map.count(pwlc_global_dof) > 0)
            {
              const auto pwld_dof = pwld_global_to_local_map.at(pwld_global_dof);
              const auto pwlc_dof = pwlc_global_to_local_map.at(pwlc_global_dof);

              pwlc_wtd_avg[pwlc_dof] += ghosted_pwld_vector[pwld_dof] * vol;
              pwlc_vol_sum[pwlc_dof] += vol;
            }
          }
  }

  // Compute nodal averages
  for (int k = 0; k < num_local_pwlc_dofs; ++k)
    pwlc_wtd_avg[k] /= pwlc_vol_sum[k];

  // Project the PWLC vector back to PWLD
  std::vector<double> output = pwld_vector;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    for (int i = 0; i < cell_mapping.GetNumNodes(); ++i)
      for (int u = 0; u < uk_man.unknowns.size(); ++u)
        for (int c = 0; c < uk_man.unknowns[u].num_components; ++c)
        {
          const auto pwld_dof = pwld.MapDOFLocal(cell, i, uk_man, u, c);
          const auto pwlc_dof = pwlc.MapDOFLocal(cell, i, uk_man, u, c);
          output[pwld_dof] = pwlc_wtd_avg[pwlc_dof];
        }
  }
  return output;
}

std::vector<int64_t>
PowerIterationKEigenSMM::MakePWLDGhostIndices(const SpatialDiscretization& pwld,
                                              const UnknownManager& uk_man)
{
  std::set<int64_t> ghost_ids;
  const auto& grid = pwld.GetGrid();
  for (const uint64_t ghost_id : grid.cells.GetGhostGlobalIDs())
  {
    const auto& cell = grid.cells[ghost_id];
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    for (int i = 0; i < cell_mapping.GetNumNodes(); ++i)
      for (int u = 0; u < uk_man.GetNumberOfUnknowns(); ++u)
        for (int c = 0; c < uk_man.unknowns[u].GetNumComponents(); ++c)
          ghost_ids.insert(pwld.MapDOF(cell, i, uk_man, u, c));
  }
  return {ghost_ids.begin(), ghost_ids.end()};
}

PowerIterationKEigenSMM::GhostInfo
PowerIterationKEigenSMM::MakePWLDGhostInfo(const SpatialDiscretization& pwld,
                                           const UnknownManager& uk_man)

{
  const auto num_local_dofs = pwld.GetNumLocalDOFs(uk_man);
  const auto num_global_dofs = pwld.GetNumGlobalDOFs(uk_man);

  // Build list of global IDs
  auto ghost_ids = MakePWLDGhostIndices(pwld, uk_man);

  // Create the ghost communicator
  auto vec_ghost_comm = std::make_shared<VectorGhostCommunicator>(
    num_local_dofs, num_global_dofs, ghost_ids, opensn::mpi_comm);

  // Create the mapping
  int64_t k = 0;
  std::map<int64_t, int64_t> ghost_global_to_local_map;
  for (const int64_t ghost_id : ghost_ids)
    ghost_global_to_local_map[ghost_id] = static_cast<int64_t>(num_local_dofs + k++);

  return {vec_ghost_comm, ghost_global_to_local_map};
}

int
PowerIterationKEigenSMM::MapAssociatedFaceNode(const Vector3& node,
                                               const std::vector<Vector3>& nbr_nodes,
                                               const double epsilon)
{
  for (int j = 0; j < nbr_nodes.size(); ++j)
    if ((node - nbr_nodes[j]).NormSquare() < epsilon)
      return j;

  throw std::logic_error("Associated neighbor node not found.");
}

} // namespace opensn
