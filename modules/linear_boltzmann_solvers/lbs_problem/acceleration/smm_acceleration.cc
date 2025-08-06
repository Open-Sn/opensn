// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT
#include "framework/object_factory.h"
#include "framework/math/parallel_vector/ghosted_parallel_stl_vector.h"
#include "framework/math/parallel_vector/parallel_stl_vector.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/smm_acceleration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/diffusion_pwlc_solver.h"
#include "modules/linear_boltzmann_solvers/solvers/pi_keigen_solver.h"

namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, SMMAcceleration);

InputParameters
SMMAcceleration::GetInputParameters()
{
  auto params = LBSKEigenAcceleration::GetInputParameters();

  params.AddOptionalParameter("sdm", "pwlc", "The spatial discretization method");

  params.ConstrainParameterRange("sdm", AllowableRangeList::New({"pwlc", "pwld"}));

  params.ChangeExistingParamToOptional("name", "SMMAcceleration");

  return params;
}

std::shared_ptr<SMMAcceleration>
SMMAcceleration::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SMMAcceleration>("lbs::SMMAcceleration", params);
}

SMMAcceleration::SMMAcceleration(const InputParameters& params)
  : LBSKEigenAcceleration(params),
    sdm_(params.GetParamValue<std::string>("sdm")),
    psi_new_local_(lbs_problem_.GetPsiNewLocal()),
    dimension_(0)
{
  lbs_problem_.GetOptions().save_angular_flux = true;
}

void
SMMAcceleration::Initialize()
{
  const auto& sdm = lbs_problem_.GetSpatialDiscretization();
  const auto num_groups = lbs_problem_.GetNumGroups();
  const auto num_gs_groups = front_gs_.groups.size();

  // Specialized SMM data
  dimension_ = lbs_problem_.GetGrid()->GetDimension();
  ComputeAuxiliaryUnitCellMatrices();
  ComputeBoundaryFactors();

  // Create PWLC structure, if needed
  if (sdm_ == "pwlc")
    InitializeLinearContinuous();

  // Create tensor structure
  for (int g = 0; g < num_groups; ++g)
    tensor_uk_man_.AddUnknown(UnknownType::VECTOR_N, dimension_ * dimension_);

  const auto local_size = sdm.GetNumLocalDOFs(tensor_uk_man_);
  const auto global_size = sdm.GetNumGlobalDOFs(tensor_uk_man_);
  const auto ghost_ids = MakePWLDGhostIndices(sdm, tensor_uk_man_);
  tensors_ = std::make_unique<GhostedParallelSTLVector>(
    local_size, global_size, ghost_ids, opensn::mpi_comm);

  // Create diffusion solver
  UnknownManager diff_uk_man;
  diff_uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

  for (const auto& [bid, bc] : lbs_problem_.GetBoundaryPreferences())
    if ((bc.type == LBSBoundaryType::ISOTROPIC) or (bc.type == LBSBoundaryType::ARBITRARY))
      throw std::logic_error("Only vacuum and reflective boundaries are valid for "
                             "k-eigenvalue problems.");

  const auto bcs = TranslateBCs(lbs_problem_.GetSweepBoundaries(), false);

  // Create the diffusion materials
  const auto xs_map = PackGroupsetXS(
    lbs_problem_.GetMatID2XSMap(), front_gs_.groups.front().id, front_gs_.groups.back().id);

  // Create the appropriate solver
  log.Log() << "Creating diffusion solver";
  if (sdm_ == "pwld")
    diffusion_solver_ = std::make_shared<DiffusionMIPSolver>(std::string(GetName() + "_SMM"),
                                                             sdm,
                                                             diff_uk_man,
                                                             bcs,
                                                             xs_map,
                                                             lbs_problem_.GetUnitCellMatrices(),
                                                             true,
                                                             verbose_);
  else
    diffusion_solver_ = std::make_shared<DiffusionPWLCSolver>(std::string(GetName() + "_SMM"),
                                                              *pwlc_ptr_,
                                                              diff_uk_man,
                                                              bcs,
                                                              xs_map,
                                                              lbs_problem_.GetUnitCellMatrices(),
                                                              true,
                                                              verbose_);

  diffusion_solver_->options.residual_tolerance = l_abs_tol_;
  diffusion_solver_->options.max_iters = max_iters_;
  diffusion_solver_->options.verbose = verbose_;
  diffusion_solver_->options.additional_options_string = petsc_options_;
  log.Log() << "Done creating diffusion solver.";

  // Initialize the solver
  log.Log() << "Initializing diffusion solver.";
  diffusion_solver_->Initialize();
  mpi_comm.barrier();
  log.Log() << "Done initializing diffusion solver.";

  // Assemble the system
  log.Log() << "Assembling diffusion system.";
  std::vector<double> tmp;
  if (sdm_ == "pwld")
    tmp.assign(sdm.GetNumLocalDOFs(diff_uk_man), 0.0);
  else
    tmp.assign(pwlc_ptr_->GetNumLocalAndGhostDOFs(diff_uk_man), 0.0);
  AssembleDiffusionBCs();
  diffusion_solver_->AssembleAand_b(tmp);
  log.Log() << "Done assembling diffusion system.";
}

void
SMMAcceleration::PreExecute()
{
  phi_ell_ = phi_old_local_;
}

void
SMMAcceleration::PrePowerIteration()
{
}

double
SMMAcceleration::PostPowerIteration()
{
  CopyOnlyPhi0(phi_new_local_, phi0_);
  phi0_old_ = phi0_;
  phi0_m_ = phi0_;

  // Update second-moment method data
  ComputeClosures(psi_new_local_);
  const auto correction = ComputeSourceCorrection();

  // Start diffusion power iterations
  double lambda = solver_->GetEigenvalue();
  double lambda_old = solver_->GetEigenvalue();
  double F0_old = lbs_problem_.ComputeFissionProduction(phi_new_local_);

  for (int m = 0; m < pi_max_its_; ++m)
  {
    // Define the nodal diffusion 1/lambda fission source.
    // This returns a local PWLD vector which must still be integrated
    // against the diffusion test functions.
    SetNodalDiffusionFissionSource(phi0_m_, Sf_);
    Scale(Sf_, 1.0 / lambda);

    // Solve the diffusion inners. This seems quite wasteful, so
    // only a single iteration is taken.
    for (int i = 0; i < 1; ++i)
    {
      // Define the nodal diffusion scattering source.
      // This returns a local PWLD vector which must still be integrated
      // against the diffusion test functions.
      SetNodalDiffusionScatterSource(phi0_old_, Ss_);

      // // Assemble the normal terms of the RHS of the diffusion solver,
      // // add in the second moment correction, and solve the system
      auto b = AssembleDiffusionRHS(Sf_ + Ss_);
      diffusion_solver_->AddToRHS(b + correction);
      diffusion_solver_->Solve(phi0_, true);

      // Bump the new solution to old
      phi0_old_ = phi0_;
    }

    // Compute a new diffusion eigenvalue. This requires mapping the
    // diffusion solution back to the transport discretization.
    std::vector<double> phi;
    ProjectBackPhi0(phi0_, phi);
    const double F0 = lbs_problem_.ComputeFissionProduction(phi);
    lambda = F0 / F0_old * lambda_old;

    // Check for convergence
    const double lambda_change = std::fabs(lambda / lambda_old - 1.0);
    const auto converged = lambda_change < pi_k_tol_;
    if (verbose_)
      log.Log() << "SMM PI iteration " << m << "  lambda " << lambda << "  change " << lambda_change
                << (converged ? "  CONVERGED" : "");

    // Bump variables for next iteration
    lambda_old = lambda;
    F0_old = F0;
    phi0_m_ = phi0_;

    if (converged)
      break;
  }

  ProjectBackPhi0(phi0_, phi_new_local_);

  const double production = lbs_problem_.ComputeFissionProduction(phi_new_local_);
  LBSVecOps::ScalePhiVector(lbs_problem_, PhiSTLOption::PHI_NEW, lambda / production);

  LBSVecOps::GSScopedCopyPrimarySTLvectors(lbs_problem_, front_gs_, phi_new_local_, phi_old_local_);
  LBSVecOps::GSScopedCopyPrimarySTLvectors(lbs_problem_, front_gs_, phi_new_local_, phi_ell_);

  return lambda;
}

void
SMMAcceleration::ComputeAuxiliaryUnitCellMatrices()
{
  const auto& discretization = lbs_problem_.GetSpatialDiscretization();

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
  const auto geom_type = lbs_problem_.GetOptions().geometry_type;
  if (geom_type == GeometryType::ONED_SPHERICAL)
    swf = std::make_shared<SphericalWeightFunction>();
  else if (geom_type == GeometryType::TWOD_CYLINDRICAL)
    swf = std::make_shared<CylindricalWeightFunction>();

  // Compute integrals
  const auto num_local_cells = lbs_problem_.GetGrid()->local_cells.size();
  K_tensor_matrices_.resize(num_local_cells);
  for (const auto& cell : lbs_problem_.GetGrid()->local_cells)
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
  if (verbose_)
    log.Log() << "Second moment method cell matrices computed.";
}

void
SMMAcceleration::ComputeBoundaryFactors()
{
  const auto& grid = lbs_problem_.GetGrid();
  const auto& pwld = lbs_problem_.GetSpatialDiscretization();
  const auto num_groupsets = lbs_problem_.GetGroupsets().size();

  // Loop over groupsets
  int gs = 0;
  for (const auto& groupset : lbs_problem_.GetGroupsets())
  {
    const auto& quad = groupset.quadrature;
    const auto num_gs_dirs = quad->omegas.size();
    const auto wt_sum = std::accumulate(quad->weights.begin(), quad->weights.end(), 0.0);

    // Loop over cells
    for (const auto& cell : grid->local_cells)
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
SMMAcceleration::AssembleDiffusionBCs() const
{
  const auto& grid = lbs_problem_.GetGrid();
  const auto& pwld = lbs_problem_.GetSpatialDiscretization();
  const auto& unit_cell_matrices = lbs_problem_.GetUnitCellMatrices();

  const auto& diff_sd = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();

  const auto num_gs_groups = front_gs_.groups.size();

  // Loop over cells
  std::vector<int64_t> rows;
  std::vector<int64_t> cols;
  std::vector<double> vals;
  for (const auto& cell : grid->local_cells)
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
        if (diffusion_solver_->GetBCS().count(face.neighbor_id))
          bc = diffusion_solver_->GetBCS().at(face.neighbor_id);

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

void
SMMAcceleration::ComputeClosures(const std::vector<std::vector<double>>& psi)
{
  const auto& grid = lbs_problem_.GetGrid();
  const auto& pwld = lbs_problem_.GetSpatialDiscretization();
  const auto& transport_views = lbs_problem_.GetCellTransportViews();

  // Create a local tensor vector, set it to zero
  auto local_tensors = tensors_->MakeLocalVector();
  local_tensors.assign(local_tensors.size(), 0.0);
  std::cerr << psi.size() << std::endl;

  // Loop over groupsets
  int gs = 0;
  for (const auto& groupset : lbs_problem_.GetGroupsets())
  {
    const auto& psi_uk_man = groupset.psi_uk_man_;
    const auto& quad = groupset.quadrature;
    const auto num_gs_dirs = quad->omegas.size();

    const auto first_grp = groupset.groups.front().id;
    const auto num_gs_groups = groupset.groups.size();

    // Loop over cells
    for (const auto& cell : grid->local_cells)
    {
      const auto& transport_view = transport_views[cell.local_id];
      const auto& cell_mapping = pwld.GetCellMapping(cell);

      // Compute node-wise, groupset wise tensors
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
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
            betas_[imap].assign(lbs_problem_.GetNumGroups(), 0.0);
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
SMMAcceleration::ComputeSourceCorrection() const
{
  const auto& grid = lbs_problem_.GetGrid();
  const auto& pwld = lbs_problem_.GetSpatialDiscretization();
  const auto& matid_to_xs_map = lbs_problem_.GetMatID2XSMap();
  const auto& unit_cell_matrices = lbs_problem_.GetUnitCellMatrices();

  const auto& diff_sd = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  auto tensors = tensors_->MakeGhostedLocalVector();

  // Create the output vector
  ParallelSTLVector output(
    diff_sd.GetNumLocalDOFs(diff_uk_man), diff_sd.GetNumGlobalDOFs(diff_uk_man), mpi_comm);

  // Build the source
  for (const auto& cell : grid->local_cells)
  {
    const auto& rho = lbs_problem_.GetDensitiesLocal()[cell.local_id];
    const auto& cell_mapping = pwld.GetCellMapping(cell);
    const auto nodes = cell_mapping.GetNodeLocations();
    const auto num_cell_nodes = cell_mapping.GetNumNodes();
    const auto num_cell_faces = cell.faces.size();

    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& K = K_tensor_matrices_[cell.local_id];

    const auto& xs = matid_to_xs_map.at(cell.block_id);
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
        const auto& nbr_cell = grid->cells[face.neighbor_id];
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
                grid->IsCellLocal(face.neighbor_id)
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

          if (sdm_ == "pwld")
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
        if (diffusion_solver_->GetBCS().count(face.neighbor_id) == 1)
          bc = diffusion_solver_->GetBCS().at(face.neighbor_id);

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
SMMAcceleration::SetNodalDiffusionFissionSource(const std::vector<double>& phi0,
                                                std::vector<double>& out) const
{
  const auto& diff_sd = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();
  const auto num_local_diff_dofs =
    pwlc_ptr_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man) : diff_sd.GetNumLocalDOFs(diff_uk_man);

  if (phi0.size() != num_local_diff_dofs)
    throw std::invalid_argument("Vector size mismatch. The flux used to construct the "
                                "diffusion fission source must have the same size as the "
                                "local diffusion solver.");

  std::vector<double> phi;
  ProjectBackPhi0(phi0, phi);

  out.resize(phi.size());
  std::fill(out.begin(), out.end(), 0.0);

  auto fn = lbs_problem_.GetActiveSetSourceFunction();
  fn(front_gs_, out, phi, APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
}

void
SMMAcceleration::SetNodalDiffusionScatterSource(const std::vector<double>& phi0,
                                                std::vector<double>& out) const
{
  const auto& diff_sd = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();
  const auto num_local_diff_dofs =
    pwlc_ptr_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man) : diff_sd.GetNumLocalDOFs(diff_uk_man);

  if (phi0.size() != num_local_diff_dofs)
    throw std::invalid_argument("Vector size mismatch. The flux used to construct the "
                                "diffusion fission source must have the same size as the "
                                "local diffusion solver.");

  std::vector<double> phi;
  ProjectBackPhi0(phi0, phi);

  out.resize(phi.size());
  std::fill(out.begin(), out.end(), 0.0);

  auto fn = lbs_problem_.GetActiveSetSourceFunction();
  fn(front_gs_,
     out,
     phi,
     APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES | SUPPRESS_WG_SCATTER);
}

std::vector<double>
SMMAcceleration::AssembleDiffusionRHS(const std::vector<double>& q0) const
{
  const auto& grid = lbs_problem_.GetGrid();
  const auto& pwld = lbs_problem_.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem_.GetUnknownManager();
  const auto& unit_cell_matrices = lbs_problem_.GetUnitCellMatrices();

  const auto& diff_sd = diffusion_solver_->GetSpatialDiscretization();
  const auto& diff_uk_man = diffusion_solver_->GetUnknownStructure();

  const auto first_grp = front_gs_.groups.front().id;
  const auto num_gs_groups = front_gs_.groups.size();

  // Clear the diffusion RHS
  const auto num_local_diff_dofs =
    pwlc_ptr_ ? diff_sd.GetNumLocalAndGhostDOFs(diff_uk_man) : diff_sd.GetNumLocalDOFs(diff_uk_man);
  std::vector<double> dummy(num_local_diff_dofs, 0.0);
  diffusion_solver_->Assemble_b(dummy);

  // Create the output vector
  ParallelSTLVector rhs(
    diff_sd.GetNumLocalDOFs(diff_uk_man), diff_sd.GetNumGlobalDOFs(diff_uk_man), mpi_comm);

  // Test the nodal source against diffusion test functions
  for (const auto& cell : grid->local_cells)
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

int
SMMAcceleration::MapAssociatedFaceNode(const Vector3& node,
                                       const std::vector<Vector3>& nbr_nodes,
                                       const double epsilon)
{
  for (int j = 0; j < nbr_nodes.size(); ++j)
    if ((node - nbr_nodes[j]).NormSquare() < epsilon)
      return j;

  throw std::logic_error("Associated neighbor node not found.");
}

} // namespace opensn
