#include "modules/mg_diffusion/mg_diffusion_solver.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
#include "framework/physics/physics_material/physics_material.h"
#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"
#include "framework/physics/physics_material/material_property_isotropic_mg_src.h"
#include "framework/physics/field_function/field_function_grid_based.h"
#include <algorithm>
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/math/spatial_discretization/finite_element/quadrature_point_data.h"
#include "modules/mg_diffusion/mg_diffusion_bndry.h"
#include "modules/mg_diffusion/tools.h"
#include <iomanip>

namespace opensn
{
namespace mg_diffusion
{

Solver::Solver(const std::string& in_solver_name)
  : opensn::Solver(in_solver_name,
                   {{"max_inner_iters", int64_t(500)},
                    {"residual_tolerance", 1.0e-2},
                    {"verbose_level", int64_t(0)},
                    {"thermal_flux_tolerance", 1.0e-2},
                    {"max_thermal_iters", int64_t(500)},
                    {"do_two_grid", false}})
{
}

Solver::~Solver()
{
  for (uint g = 0; g < num_groups_; ++g)
  {
    VecDestroy(&x_[g]);
    VecDestroy(&bext_[g]);
    MatDestroy(&A_[g]);
  }
  VecDestroy(&b_);

  if (last_fast_group_ < num_groups_)
  {
    VecDestroy(&thermal_dphi_);
    for (uint g = last_fast_group_; g < num_groups_; ++g)
      VecDestroy(&x_old_[g]);
  }

  if (do_two_grid_)
  {
    VecDestroy(&x_[num_groups_]);
    MatDestroy(&A_[num_groups_]);
  }
}

void
Solver::Initialize()
{
  log.Log() << "\n"
            << program_timer.GetTimeString() << " " << TextName()
            << ": Initializing CFEM Multigroup Diffusion solver ";

  // Get grid
  grid_ptr_ = GetCurrentHandler().GetGrid();
  const auto& grid = *grid_ptr_;
  if (grid_ptr_ == nullptr)
    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + " No grid defined.");

  log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Add unique material ids
  std::set<int> unique_material_ids;
  int invalid_mat_cell_count = 0;
  for (auto& cell : grid.local_cells)
  {
    unique_material_ids.insert(cell.material_id_);
    if (cell.material_id_ < 0) ++invalid_mat_cell_count;
  }
  const auto& ghost_cell_ids = grid.cells.GetGhostGlobalIDs();
  for (uint64_t cell_id : ghost_cell_ids)
  {
    const auto& cell = grid.cells[cell_id];
    unique_material_ids.insert(cell.material_id_);
    if (cell.material_id_ < 0) ++invalid_mat_cell_count;
  }

  if (invalid_mat_cell_count > 0)
  {
    log.LogAllWarning() << "Number of invalid material cells: " << invalid_mat_cell_count;
  }

  // Initialize materials
  mg_diffusion::Solver::Initialize_Materials(unique_material_ids);

  // BIDs
  auto globl_unique_bndry_ids = grid.GetDomainUniqueBoundaryIDs();
  mg_diffusion::Solver::Set_BCs(globl_unique_bndry_ids);

  // Make SDM
  sdm_ptr_ = PieceWiseLinearContinuous::New(*grid_ptr_);
  const auto& sdm = *sdm_ptr_;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;
  num_local_dofs_ = sdm.GetNumLocalDOFs(OneDofPerNode);
  num_globl_dofs_ = sdm.GetNumGlobalDOFs(OneDofPerNode);

  log.Log() << "Num local DOFs: " << num_local_dofs_;
  log.Log() << "Num globl DOFs: " << num_globl_dofs_;

  // Initializes Mats and Vecs
  const auto n = static_cast<int64_t>(num_local_dofs_);
  const auto N = static_cast<int64_t>(num_globl_dofs_);

  std::vector<int64_t> nodal_nnz_in_diag;
  std::vector<int64_t> nodal_nnz_off_diag;
  sdm.BuildSparsityPattern(nodal_nnz_in_diag, nodal_nnz_off_diag, OneDofPerNode);

  unsigned int i_two_grid = do_two_grid_ ? 1 : 0;
  //  std::cout << "i_two_grid = " << i_two_grid << std::endl;

  A_.resize(num_groups_ + i_two_grid, nullptr);
  x_.resize(num_groups_ + i_two_grid, nullptr);
  bext_.resize(num_groups_, nullptr);

  auto ghost_dof_indices = sdm.GetGhostDOFIndices(OneDofPerNode);

  for (uint g = 0; g < num_groups_; ++g)
  {
    // x[g] = CreateVector(n,N);
    x_[g] = CreateVectorWithGhosts(
      n, N, static_cast<int64_t>(ghost_dof_indices.size()), ghost_dof_indices);
    VecSet(x_[g], 0.0);
    bext_[g] = CreateVector(n, N);

    A_[g] = CreateSquareMatrix(n, N);
    InitMatrixSparsity(A_[g], nodal_nnz_in_diag, nodal_nnz_off_diag);
  }
  // initialize b
  VecDuplicate(bext_.front(), &b_);
  // initialize old flux for thermal groups
  if (last_fast_group_ < num_groups_)
  {
    VecDuplicate(x_.front(), &thermal_dphi_);
    x_old_.resize(num_groups_, nullptr);
    for (uint g = 0; g < num_groups_; ++g)
    {
      VecDuplicate(x_.front(), &x_old_[g]); // jcr is x_old like bext or like x?
      VecSet(x_old_[g], 0.0);
    }
  }
  // add two-grid mat and vec, if needed
  if (do_two_grid_)
  {
    A_[num_groups_] = CreateSquareMatrix(n, N);
    InitMatrixSparsity(A_[num_groups_], nodal_nnz_in_diag, nodal_nnz_off_diag);
    VecDuplicate(x_.front(), &x_[num_groups_]); // jcr needed?
    //    x[num_groups] = CreateVectorWithGhosts(n,N,
    //                                                        static_cast<int64_t>(ghost_dof_indices.size()),
    //                                                        ghost_dof_indices);
  }

  if (do_two_grid_) mg_diffusion::Solver::Compute_TwoGrid_VolumeFractions();

  // Create Mats and ExtVecs
  mg_diffusion::Solver::Assemble_A_bext();

  // Volume fraction for two-grid update

  // Field Function
  if (field_functions_.empty())
  {
    for (uint g = 0; g < mg_diffusion::Solver::num_groups_; ++g)
    {
      std::string solver_name;
      if (not TextName().empty()) solver_name = TextName() + "-";

      char buff[100];
      int dummy = snprintf(buff, 4, "%03d", g);

      std::string text_name = solver_name + "Flux_g" + std::string(buff);

      auto initial_field_function =
        std::make_shared<FieldFunctionGridBased>(text_name, sdm_ptr_, Unknown(UnknownType::SCALAR));

      field_functions_.push_back(initial_field_function);
      field_function_stack.push_back(initial_field_function);
    } // for g
  }   // if not ff set
}

void
Solver::Initialize_Materials(std::set<int>& material_ids)
{
  log.Log0Verbose1() << "Initializing Materials";

  std::stringstream materials_list;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process materials found
  const size_t num_physics_mats = material_stack.size();
  bool first_material_read = true;

  for (const int& mat_id : material_ids)
  {
    auto current_material = GetStackItemPtr(material_stack, mat_id, __FUNCTION__);
    materials_list << "Material id " << mat_id;

    // Check valid ids
    if (mat_id < 0)
    {
      throw std::logic_error("MG-diff-InitMaterials: Cells encountered with no assigned material.");
    }
    if (static_cast<size_t>(mat_id) >= num_physics_mats)
    {
      throw std::logic_error("MG-diff-InitMaterials: Cells encountered with "
                             "material id that matches no material in physics material library.");
    }

    // Extract properties
    using MatProperty = PropertyType;
    bool found_transport_xs = false;
    for (const auto& property : current_material->properties_)
    {
      if (property->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs = std::static_pointer_cast<MultiGroupXS>(property);
        matid_to_xs_map[mat_id] = transp_xs;
        found_transport_xs = true;
        if (first_material_read) num_groups_ = transp_xs->NumGroups();

      } // transport xs
      if (property->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source = std::static_pointer_cast<IsotropicMultiGrpSource>(property);

        if (mg_source->source_value_g_.size() < num_groups_)
        {
          log.LogAllWarning() << "MG-Diff-InitMaterials: Isotropic Multigroup source specified "
                                 "in "
                              << "material \"" << current_material->name_ << "\" has fewer "
                              << "energy groups than called for in the simulation. "
                              << "Source will be ignored.";
        }
        else { matid_to_src_map[mat_id] = mg_source; }
      } // P0 source
    }   // for property

    // Check valid property
    if (!found_transport_xs)
    {
      log.LogAllError() << "MG-Diff-InitMaterials: Found no transport cross-section property "
                           "for "
                        << "material \"" << current_material->name_ << "\".";
      Exit(EXIT_FAILURE);
    }
    // Check number of groups legal
    if (matid_to_xs_map[mat_id]->NumGroups() != num_groups_)
    {
      log.LogAllError() << "MG-Diff-InitMaterials: Found material \"" << current_material->name_
                        << "\" has " << matid_to_xs_map[mat_id]->NumGroups() << " groups and "
                        << "the simulation has " << num_groups_ << " groups. The material "
                        << "must have the same number of groups.";
      Exit(EXIT_FAILURE);
    }

    // Check number of moments
    if (matid_to_xs_map[mat_id]->ScatteringOrder() > 1)
    {
      log.Log0Warning() << "MG-Diff-InitMaterials: Found material \"" << current_material->name_
                        << "\" has a scattering order of "
                        << matid_to_xs_map[mat_id]->ScatteringOrder() << " and"
                        << " the simulation has a scattering order of One (MG-Diff)"
                        << " The higher moments will therefore not be used.";
    }

    materials_list << " number of moments " << matid_to_xs_map[mat_id]->ScatteringOrder() + 1
                   << "\n";

    first_material_read = false;
  } // for material id

  log.Log() << "Materials Initialized:\n" << materials_list.str() << "\n";

  opensn::mpi_comm.barrier();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute last fast group
  // initialize last fast group
  log.Log() << "Computing last fast group.";
  unsigned int lfg = num_groups_;

  if (num_groups_ > 1)
  {
    // loop over all materials
    for (const auto& mat_id_xs : matid_to_xs_map)
    {
      // get the P0 transfer matrix
      const auto& S = mat_id_xs.second->TransferMatrix(0);
      // loop over all row of the transfer matrix
      const int G = static_cast<int>(num_groups_);
      for (int g = G - 1; g >= 0; --g)
      {
        for (const auto& [row_g, gp, sigma_sm] : S.Row(g))
        {
          if ((std::fabs(sigma_sm) > 1e-10) && (gp > row_g))
            lfg = std::min(lfg, static_cast<unsigned int>(row_g));
        }
      }
    } // loop on materials
  }   // if num_groups>1

  last_fast_group_ = lfg;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute two-grid params
  do_two_grid_ = basic_options_("do_two_grid").BoolValue();
  if ((lfg == num_groups_) && do_two_grid_)
  {
    log.Log0Error() << "Two-grid is not possible with no upscattering.";
    do_two_grid_ = false;
    Exit(EXIT_FAILURE);
  }
  if (do_two_grid_)
  {
    log.Log() << "Compute_TwoGrid_Params";
    Compute_TwoGrid_Params();
  }
}

void
Solver::Compute_TwoGrid_Params()
{
  // loop over all materials
  for (const auto& mat_id_xs : matid_to_xs_map)
  {

    // get the P0 transfer matrix and total XS
    const auto& isotropic_transfer_matrix = mat_id_xs.second->TransferMatrix(0);
    const auto& sigma_t = mat_id_xs.second->SigmaTotal();
    const auto& diffusion_coeff = mat_id_xs.second->DiffusionCoefficient();

    // put P0 transfer matrix in nicer form
    MatDbl S(num_groups_, VecDbl(num_groups_, 0.0));
    for (unsigned int g = 0; g < num_groups_; ++g)
      for (const auto& [row_g, gprime, sigma] : isotropic_transfer_matrix.Row(g))
        S[g][gprime] = sigma;

    // (L+D) e_new = -U e_old
    // original matrix = diag(total) - scattering
    // so L+D = diag(removal) - tril(scattering)
    // and U = -triu(scattering)
    MatDbl A(num_groups_, VecDbl(num_groups_, 0.0));
    MatDbl B(num_groups_, VecDbl(num_groups_, 0.0));
    for (unsigned int g = 0; g < num_groups_; ++g)
    {
      A[g][g] = sigma_t[g] - S[g][g];
      for (unsigned int gp = 0; gp < g; ++gp)
        A[g][gp] = -S[g][gp];
      for (unsigned int gp = g + 1; gp < num_groups_; ++gp)
        B[g][gp] = S[g][gp];
    }
    MatDbl Ainv = Inverse(A);
    // finally, obtain the iteration matrix
    MatDbl C_ = MatMul(Ainv, B);
    // Perform power iteration
    VecDbl E(num_groups_, 1.0);
    double rho = PowerIteration(C_, E, 10000, 1.0e-12);

    // Compute two-grid diffusion quantities
    // normalize spectrum
    std::vector<double> spectrum(num_groups_, 1.0);
    double sum = 0.0;
    for (unsigned int g = 0; g < num_groups_; ++g)
      sum += std::fabs(E[g]);
    for (unsigned int g = 0; g < num_groups_; ++g)
      spectrum[g] = std::fabs(E[g]) / sum;
    // D ave and Sigma_a ave
    double collapsed_D = 0.0;
    double collapsed_sig_a = 0.0;
    for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
    {
      collapsed_D += diffusion_coeff[g] * spectrum[g];
      collapsed_sig_a += sigma_t[g] * spectrum[g];
      for (unsigned int gp = last_fast_group_; gp < num_groups_; ++gp)
        collapsed_sig_a -= S[g][gp] * spectrum[gp];
    }
    // Verbose output the spectrum
    log.Log0Verbose1() << "Fundamental eigen-value: " << rho;
    std::stringstream outstr;
    for (auto& xi : spectrum)
      outstr << xi << '\n';
    log.Log0Verbose1() << outstr.str(); // jcr verbose1

    //    std::stringstream outstr2;
    //    for (auto &xi: diffusion_coeff)
    //      outstr2 << xi << '\n';
    //    chi::log.Log0Verbose0() << outstr2.str();  // jcr verbose1
    //
    //    std::cout << "collapsed = " << collapsed_sig_a
    //    <<", "<< collapsed_D << std::endl;
    //    chi::Exit(12345);

    const auto mat_id = mat_id_xs.first;
    map_mat_id_2_tginfo.insert(
      std::make_pair(mat_id, TwoGridCollapsedInfo{collapsed_D, collapsed_sig_a, spectrum}));

  } // end loop over materials
}

void
Solver::Compute_TwoGrid_VolumeFractions()
{
  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  const size_t ncells = grid.local_cells.size();
  VF_.resize(ncells);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    VF_[counter].resize(num_nodes, 0.0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      double vol_frac_shape_i = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        vol_frac_shape_i += qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      vol_frac_shape_i /= cell_mapping.CellVolume();
      VF_[counter][i] = vol_frac_shape_i;
    } // for i

    counter++;
  } // for cell
}

void
Solver::Set_BCs(const std::vector<uint64_t>& globl_unique_bndry_ids)
{
  log.Log0Verbose1() << "Setting Boundary Conditions";

  uint64_t max_boundary_id = 0;
  for (const auto& id : globl_unique_bndry_ids)
    max_boundary_id = std::max(id, max_boundary_id);

  log.Log() << "Max boundary id identified: " << max_boundary_id;

  for (int bndry = 0; bndry < (max_boundary_id + 1); ++bndry)
  {
    if (boundary_preferences_.find(bndry) != boundary_preferences_.end())
    {
      BoundaryInfo bndry_info = boundary_preferences_.at(bndry);
      auto& bndry_vals = bndry_info.second;

      switch (bndry_info.first)
      {
        case BoundaryType::Reflecting: // ------------- REFLECTING
        {
          boundaries_.push_back({BoundaryType::Reflecting});
          log.Log() << "Boundary " << bndry << " set to reflecting.";
          break;
        }
        case BoundaryType::Robin: // ------------- ROBIN
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Robin needs 3 values in bndry vals.");
          boundaries_.push_back(Boundary{BoundaryType::Robin, bndry_vals});
          log.Log() << "Boundary " << bndry << " set to robin.";
          break;
        }
        case BoundaryType::Vacuum: // ------------- VACUUM
        {
          auto ng = mg_diffusion::Solver::num_groups_;
          std::vector<double> a_values(ng, 0.25);
          std::vector<double> b_values(ng, 0.5);
          std::vector<double> f_values(ng, 0.0);
          boundaries_.push_back(Boundary{BoundaryType::Robin, {a_values, b_values, f_values}});
          log.Log() << "Boundary " << bndry << " set to vacuum.";
          break;
        }
        case BoundaryType::Neumann: // ------------- NEUMANN
        {
          if (bndry_vals.size() != 3)
            throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                   " Neumann needs 3 values in bndry vals.");
          boundaries_.push_back(Boundary{BoundaryType::Robin, bndry_vals});
          log.Log() << "Boundary " << bndry << " set to neumann.";
          break;
        }
      } // switch boundary type
    }
    else
    {
      auto ng = mg_diffusion::Solver::num_groups_;
      std::vector<double> a_values(ng, 0.25);
      std::vector<double> b_values(ng, 0.5);
      std::vector<double> f_values(ng, 0.0);
      boundaries_.push_back({BoundaryType::Robin, {a_values, b_values, f_values}});
      log.Log0Verbose1() << "No boundary preference found for boundary index " << bndry
                         << "Vacuum boundary added as default.";
    }
  } // for bndry
}

void
Solver::Assemble_A_bext()
{
  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  // Assemble the system
  unsigned int i_two_grid = do_two_grid_ ? 1 : 0;

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& xs = matid_to_xs_map.at(cell.material_id_);
    const auto& D = xs->DiffusionCoefficient();
    const auto& sigma_r = xs->SigmaRemoval();

    const auto& qext = matid_to_src_map.at(cell.material_id_);
    double collapsed_D = 0., collapsed_sig_a = 0.;
    if (do_two_grid_)
    {
      const auto& xstg = map_mat_id_2_tginfo.at(cell.material_id_);
      collapsed_D = xstg.collapsed_D;
      collapsed_sig_a = xstg.collapsed_sig_a;
    }

    std::vector<VecDbl> rhs_cell;
    rhs_cell.resize(num_groups_);
    for (uint g = 0; g < num_groups_; ++g)
      rhs_cell[g].resize(num_nodes, 0.0);

    std::vector<MatDbl> Acell;
    Acell.resize(num_groups_ + i_two_grid);
    for (uint g = 0; g < num_groups_ + i_two_grid; ++g)
      Acell[g].resize(num_nodes, VecDbl(num_nodes, 0.0));

    for (size_t i = 0; i < num_nodes; ++i)
    {
      for (size_t j = 0; j < num_nodes; ++j)
      {
        double entry_mij = 0.0;
        double entry_kij = 0.0;
        for (size_t qp : qp_data.QuadraturePointIndices())
        {
          entry_mij += qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp) * qp_data.JxW(qp);

          entry_kij += qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp)) * qp_data.JxW(qp);
        } // for qp
        for (uint g = 0; g < num_groups_; ++g)
          Acell[g][i][j] = entry_mij * sigma_r[g] + entry_kij * D[g];

        if (do_two_grid_)
          Acell[num_groups_][i][j] = entry_mij * collapsed_sig_a + entry_kij * collapsed_D;
      } // for j
      double entry_rhsi = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        entry_rhsi += qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      for (uint g = 0; g < num_groups_; ++g)
        rhs_cell[g][i] = entry_rhsi * (qext->source_value_g_[g]);
    } // for i

    // Deal with BC (all based on variations of Robin)
    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      // not a boundary face
      if (face.has_neighbor_) continue;

      auto& bndry = boundaries_[face.neighbor_id_];

      // Robin boundary
      //   for two-grid, it is homogenous Robin
      if (bndry.type_ == BoundaryType::Robin)
      {
        const auto qp_face_data = cell_mapping.MakeSurfaceQuadraturePointData(f);
        const size_t num_face_nodes = face.vertex_ids_.size();

        auto& aval = bndry.mg_values_[0];
        auto& bval = bndry.mg_values_[1];
        auto& fval = bndry.mg_values_[2];
        if (do_two_grid_)
        {
          aval.push_back(0.25);
          bval.push_back(0.5);
          fval.push_back(0.0);
        }

        // sanity check, Assert if b=0
        for (uint g = 0; g < num_groups_ + i_two_grid; ++g)
        {
          if (std::fabs(bval[g]) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");
        }

        // true Robin when a!=0, otherwise, it is a Neumann:
        // only do this part if true Robin (i.e., a!=0)
        for (uint g = 0; g < num_groups_ + i_two_grid; ++g)
        {
          if (std::fabs(aval[g]) > 1.0e-8)
          {
            // loop over nodes of that face
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const uint i = cell_mapping.MapFaceNode(f, fi);

              double entry_rhsi = 0.0;
              for (size_t qp : qp_face_data.QuadraturePointIndices())
                entry_rhsi += qp_face_data.ShapeValue(i, qp) * qp_face_data.JxW(qp);
              if (g < num_groups_) // check due to two-grid
                rhs_cell[g][i] += fval[g] / bval[g] * entry_rhsi;

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const uint j = cell_mapping.MapFaceNode(f, fj);
                double entry_aij = 0.0;
                for (size_t qp : qp_face_data.QuadraturePointIndices())
                  entry_aij += qp_face_data.ShapeValue(i, qp) * qp_face_data.ShapeValue(j, qp) *
                               qp_face_data.JxW(qp);
                Acell[g][i][j] += aval[g] / bval[g] * entry_aij;
              } // for fj
            }   // for fi
          }     // end true Robin
        }       // for g
      }         // if Robin
    }           // for face f

    // Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); // node-mapping
    for (size_t i = 0; i < num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    // Assembly into system
    for (uint g = 0; g < num_groups_; ++g)
      for (size_t i = 0; i < num_nodes; ++i)
        VecSetValue(bext_[g], imap[i], rhs_cell[g][i], ADD_VALUES);

    for (uint g = 0; g < num_groups_ + i_two_grid; ++g)
      for (size_t i = 0; i < num_nodes; ++i)
        for (size_t j = 0; j < num_nodes; ++j)
          MatSetValue(A_[g], imap[i], imap[j], Acell[g][i][j], ADD_VALUES);

  } // for cell

  log.Log() << "Global assembly";

  for (uint g = 0; g < num_groups_; ++g)
  {
    VecAssemblyBegin(bext_[g]);
    VecAssemblyEnd(bext_[g]);
  }
  for (uint g = 0; g < num_groups_ + i_two_grid; ++g)
  {
    MatAssemblyBegin(A_[g], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_[g], MAT_FINAL_ASSEMBLY);
  }

  //  PetscViewer viewer;
  //  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A2_before_bc.m",&viewer);
  //  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  //  MatView(A[0],viewer);
  //  PetscViewerPopFormat(viewer);
  //  PetscViewerDestroy(&viewer);
  //  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"bext2_before_bc.m",&viewer);
  //  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  //  VecView(bext[0],viewer);
  //  PetscViewerPopFormat(viewer);
  //  PetscViewerDestroy(&viewer);

  log.Log() << "Done global assembly";
}

void
Solver::Execute()
{
  log.Log() << "\nExecuting CFEM Multigroup Diffusion solver";

  // Create Krylov Solver
  // setup KSP once for all
  petsc_solver_ = CreateCommonKrylovSolverSetup(A_.front(),
                                                TextName(),
                                                KSPCG,
                                                PCGAMG,
                                                basic_options_("residual_tolerance").FloatValue(),
                                                basic_options_("max_inner_iters").IntegerValue());

  KSPSetApplicationContext(petsc_solver_.ksp, (void*)&my_app_context_);
  KSPMonitorCancel(petsc_solver_.ksp);
  KSPMonitorSet(petsc_solver_.ksp, &mg_diffusion::MGKSPMonitor, nullptr, nullptr);

  int64_t iverbose = basic_options_("verbose_level").IntegerValue();
  my_app_context_.verbose = iverbose > 1 ? PETSC_TRUE : PETSC_FALSE;
  //  if (my_app_context.verbose == PETSC_TRUE)
  //    cout << "--context TRUE" << endl;
  //  if (my_app_context.verbose == PETSC_FALSE)
  //    cout << "--context FALSE" << endl;
  //  std::cout << "STOP" << std::endl; std::cin.get();

  // shortcuts
  // unsigned int lfg = mg_diffusion::Solver::last_fast_group;
  // unsigned int ng = mg_diffusion::Solver::num_groups;

  // Solve fast groups:
  for (unsigned int g = 0; g < last_fast_group_; ++g)
  {
    mg_diffusion::Solver::Assemble_RHS(g, iverbose);
    mg_diffusion::Solver::SolveOneGroupProblem(g, iverbose);
  }

  // Solve thermal groups:
  unsigned int thermal_iteration = 0;
  // max # of thermal iterations
  int64_t max_thermal_iters = basic_options_("max_thermal_iters").IntegerValue();
  // max thermal error between two successive iterates
  double thermal_tol = basic_options_("thermal_flux_tolerance").FloatValue();
  // computed error
  double thermal_error_all;
  double thermal_error_g;

  do
  {
    thermal_error_all = 0.0;
    for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
    {
      // conpute rhs src
      mg_diffusion::Solver::Assemble_RHS(g, iverbose);
      // copy solution
      VecCopy(x_[g], x_old_[g]);
      // solve group g for new solution
      mg_diffusion::Solver::SolveOneGroupProblem(g, iverbose);
      // compute L2 norm of thermal error for current g (requires one more copy)
      VecCopy(x_[g], thermal_dphi_);
      VecAXPY(thermal_dphi_, -1.0, x_old_[g]);
      VecNorm(thermal_dphi_, NORM_2, &thermal_error_g);
      thermal_error_all = std::max(thermal_error_all, thermal_error_g);
    }
    // perform two-grid
    if (do_two_grid_)
    {
      mg_diffusion::Solver::Assemble_RHS_TwoGrid(iverbose);
      mg_diffusion::Solver::SolveOneGroupProblem(num_groups_, iverbose);
      mg_diffusion::Solver::Update_Flux_With_TwoGrid(iverbose);
    }

    if (iverbose > 0)
      log.Log() << " --thermal iteration = " << std::setw(5) << std::right << thermal_iteration
                << ", Error=" << std::setw(11) << std::right << std::scientific
                << std::setprecision(7) << thermal_error_all << std::endl;

    ++thermal_iteration;
  } while ((thermal_error_all > thermal_tol) && (thermal_iteration < max_thermal_iters));

  if (iverbose > 0)
  {
    if (thermal_error_all < thermal_tol)
      std::cout << "\nThermal iterations converged for fixed-source problem" << std::endl;
    else
      std::cout << "\nThermal iterations NOT converged for fixed-source problem" << std::endl;
  }

  UpdateFieldFunctions();
  log.Log() << "Done solving multi-group diffusion";
}

void
Solver::Assemble_RHS(const unsigned int g, const int64_t verbose)
{
  if (verbose > 2) log.Log() << "\nAssemblying RHS for group " + std::to_string(g);

  // copy the external source vector for group g into b
  VecSet(b_, 0.0);
  VecCopy(bext_[g], b_);

  const auto& sdm = *mg_diffusion::Solver::sdm_ptr_;
  // compute inscattering term
  for (const auto& cell : mg_diffusion::Solver::grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& xs = matid_to_xs_map.at(cell.material_id_);
    const auto& S = xs->TransferMatrix(0);

    for (const auto& [row_g, gprime, sigma_sm] : S.Row(g))
    {
      if (gprime != g) // g and row_g are the same, maybe different int types
      {
        const double* xlocal;
        VecGetArrayRead(x_[gprime], &xlocal);

        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOF(cell, i);
          double inscatter_g = 0.0;

          for (size_t j = 0; j < num_nodes; ++j)
          {
            const int64_t jmap = sdm.MapDOFLocal(cell, j);

            // get flux at node j
            const double flxj_gp = xlocal[jmap];
            for (size_t qp : qp_data.QuadraturePointIndices())
              inscatter_g += sigma_sm * flxj_gp * qp_data.ShapeValue(i, qp) *
                             qp_data.ShapeValue(j, qp) * qp_data.JxW(qp);
          } // for j
          // add inscattering value to vector
          VecSetValue(b_, imap, inscatter_g, ADD_VALUES);
        } // for i
        VecRestoreArrayRead(x_[gprime], &xlocal);
      } // if gp!=g
    }   // for gprime
  }     // for cell

  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);
}

void
Solver::SolveOneGroupProblem(const unsigned int g, const int64_t verbose)
{
  if (verbose > 1) log.Log() << "Solving group: " << g;

  KSPSetOperators(petsc_solver_.ksp, A_[g], A_[g]);
  KSPSolve(petsc_solver_.ksp, b_, x_[g]);

  // this is required to compute the inscattering RHS correctly in parallel
  CommunicateGhostEntries(x_[g]);

  if (verbose > 1) log.Log() << "Done solving group " << g;
}

void
Solver::Assemble_RHS_TwoGrid(const int64_t verbose)
{
  if (verbose > 2) log.Log() << "\nAssemblying RHS for two-grid ";

  VecSet(b_, 0.0);

  const auto& sdm = *sdm_ptr_;
  // compute inscattering term
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumetricQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& S = matid_to_xs_map.at(cell.material_id_)->TransferMatrix(0);

    for (unsigned g = last_fast_group_; g < num_groups_; ++g)
    {
      for (const auto& [row_g, gprime, sigma_sm] : S.Row(g))
      {
        if (gprime > g) // the upper part for the residual of two-grid accel
        {
          const double* xlocal;
          const double* xlocal_old;
          VecGetArrayRead(x_[gprime], &xlocal);
          VecGetArrayRead(x_old_[gprime], &xlocal_old);

          for (size_t i = 0; i < num_nodes; ++i)
          {
            const int64_t imap = sdm.MapDOF(cell, i);
            double inscatter_g = 0.0;

            for (size_t j = 0; j < num_nodes; ++j)
            {
              const int64_t jmap = sdm.MapDOFLocal(cell, j);

              // get flux at node j
              const double delta_flxj_gp = xlocal[jmap] - xlocal_old[jmap];
              for (size_t qp : qp_data.QuadraturePointIndices())
                inscatter_g += sigma_sm * delta_flxj_gp * qp_data.ShapeValue(i, qp) *
                               qp_data.ShapeValue(j, qp) * qp_data.JxW(qp);
            } // for j
            // add inscattering value to vector
            VecSetValue(b_, imap, inscatter_g, ADD_VALUES);
          } // for i
          VecRestoreArrayRead(x_[gprime], &xlocal);
          VecRestoreArrayRead(x_old_[gprime], &xlocal_old);
        } // if gp!=g
      }   // for gprime
    }     // for g
  }       // for cell

  VecAssemblyBegin(b_);
  VecAssemblyEnd(b_);
  //  VecView(b, PETSC_VIEWER_STDERR_WORLD);
  //  chi::Exit(1234);
}

void
Solver::Update_Flux_With_TwoGrid(const int64_t verbose)
{
  if (verbose > 2) log.Log() << "\nUpdating Thermal fluxes from two-grid";

  const auto& grid = *grid_ptr_;
  const auto& sdm = *sdm_ptr_;

  // contains two_grid flux, stored in last num_groups entry
  const double* xlocal_tg;
  VecGetArrayRead(x_[num_groups_], &xlocal_tg);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& xstg = map_mat_id_2_tginfo.at(cell.material_id_);

    for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
    {
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imap = sdm.MapDOFLocal(cell, i);
        const int64_t imap_globl = sdm.MapDOF(cell, i);
        const double aux = xlocal_tg[imap] * VF_[counter][i] * xstg.spectrum[g];
        VecSetValue(x_[g], imap_globl, aux, ADD_VALUES);
      } // i
    }   // g
    counter++;
  } // for cell

  // finalize
  for (unsigned int g = last_fast_group_; g < num_groups_; ++g)
  {
    VecAssemblyBegin(x_[g]);
    VecAssemblyEnd(x_[g]);
  }
  // release two-grid flux
  VecRestoreArrayRead(x_[num_groups_], &xlocal_tg);
}

void
Solver::UpdateFieldFunctions()
{
  const auto& OneDOFPerNode = sdm_ptr_->UNITARY_UNKNOWN_MANAGER;
  for (int g = 0; g < num_groups_; ++g)
  {
    std::vector<double> data_vector;
    sdm_ptr_->LocalizePETScVector(x_[g], data_vector, OneDOFPerNode);

    auto& ff = field_functions_.at(g);
    ff->UpdateFieldVector(data_vector);
  } // for g
}

} // namespace mg_diffusion
} // namespace opensn
