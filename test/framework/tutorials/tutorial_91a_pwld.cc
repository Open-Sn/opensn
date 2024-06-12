#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "lua/framework/console/console.h"
#include "framework/math/math_range.h"
#include "framework/data_types/ndarray.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>

using namespace opensn;

namespace unit_sim_tests
{

/**PWLD Sweep. */
ParameterBlock SimTest91_PWLD(const InputParameters&);

RegisterWrapperFunctionInNamespace(unit_tests, SimTest91_PWLD, nullptr, SimTest91_PWLD);

ParameterBlock
SimTest91_PWLD(const InputParameters&)
{
  const std::string fname = "SimTest91_PWLD";

  opensn::log.Log() << "SimTest91_PWLD num_args = " << 0;

  if (opensn::mpi_comm.size() != 1)
    throw std::logic_error(fname + ": Is serial only.");

  // Get grid
  auto grid_ptr = GetCurrentMesh();
  const auto& grid = *grid_ptr;

  opensn::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make Orthogonal mapping
  const auto ijk_info = grid.GetIJKInfo();
  const auto& ijk_mapping = grid.MakeIJKToGlobalIDMapping();

  const auto Nx = static_cast<int64_t>(ijk_info[0]);
  const auto Ny = static_cast<int64_t>(ijk_info[1]);
  const auto Nz = static_cast<int64_t>(ijk_info[2]);

  auto dimension = grid.Dimension();

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

  const size_t num_local_nodes = sdm.GetNumLocalDOFs(OneDofPerNode);
  const size_t num_globl_nodes = sdm.GetNumGlobalDOFs(OneDofPerNode);

  opensn::log.Log() << "Num local nodes: " << num_local_nodes;
  opensn::log.Log() << "Num globl nodes: " << num_globl_nodes;

  // Make an angular quadrature
  std::shared_ptr<AngularQuadrature> quadrature;
  if (dimension == 1)
    quadrature = std::make_shared<AngularQuadratureProdGL>(8);
  else if (dimension == 2)
  {
    quadrature = std::make_shared<AngularQuadratureProdGLC>(8, 8);
    quadrature->OptimizeForPolarSymmetry(4.0 * M_PI);
  }
  else if (dimension == 3)
    quadrature = std::make_shared<AngularQuadratureProdGLC>(8, 8);
  else
    throw std::logic_error(fname + "Error with the dimensionality "
                                   "of the mesh.");
  opensn::log.Log() << "Quadrature created." << std::endl;

  // Set/Get params
  const size_t scat_order = 1;
  const size_t num_groups = 20;

  quadrature->BuildMomentToDiscreteOperator(scat_order, dimension);
  quadrature->BuildDiscreteToMomentOperator(scat_order, dimension);

  const auto& m2d = quadrature->GetMomentToDiscreteOperator();
  const auto& d2m = quadrature->GetDiscreteToMomentOperator();
  const auto& m_ell_em_map = quadrature->GetMomentToHarmonicsIndexMap();

  const size_t num_moments = m_ell_em_map.size();
  const size_t num_dirs = quadrature->omegas_.size();

  opensn::log.Log() << "End Set/Get params." << std::endl;
  opensn::log.Log() << "Num Moments: " << num_moments << std::endl;

  // Make Unknown Managers
  const auto VecN = UnknownType::VECTOR_N;
  using Unknown = Unknown;

  std::vector<Unknown> phi_uks(num_moments, Unknown(VecN, num_groups));
  std::vector<Unknown> psi_uks(num_dirs, Unknown(VecN, num_groups));

  const UnknownManager phi_uk_man(phi_uks);
  const UnknownManager psi_uk_man(psi_uks);

  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(phi_uk_man);
  const size_t num_local_psi_dofs = sdm.GetNumLocalDOFs(psi_uk_man);

  opensn::log.Log() << "End ukmanagers." << std::endl;

  // Make XSs
  MultiGroupXS xs;
  xs.Initialize("xs_graphite_pure.xs");

  // Initializes vectors
  std::vector<double> phi_old(num_local_phi_dofs, 0.0);
  std::vector<double> psi_old(num_local_psi_dofs, 0.0);
  auto source_moments = phi_old;
  auto phi_new = phi_old;
  auto q_source = phi_old;

  opensn::log.Log() << "End vectors." << std::endl;

  // Make material source term
  for (const auto& cell : grid.local_cells)
  {
    const auto& cc = cell.centroid_;
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    if (cc.x < 0.5 and cc.y < 0.5 and cc.z < 0.5 and cc.x > -0.5 and cc.y > -0.5 and cc.z > -0.5)
    {
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

        q_source[dof_map] = 1.0;
      }
    }
  }

  // Precompute cell matrices
  typedef std::vector<Vector3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;
  typedef std::vector<MatDbl> VecMatDbl;

  std::vector<MatVec3> cell_Gmatrices;
  std::vector<MatDbl> cell_Mmatrices;
  std::vector<VecMatDbl> cell_faceMmatrices;

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    MatVec3 IntV_shapeI_gradshapeJ(num_nodes, VecVec3(num_nodes, Vector3(0, 0, 0)));
    MatDbl IntV_shapeI_shapeJ(num_nodes, VecDbl(num_nodes, 0.0));

    for (unsigned int i = 0; i < num_nodes; ++i)
      for (unsigned int j = 0; j < num_nodes; ++j)
        for (const auto& qp : fe_vol_data.QuadraturePointIndices())
        {
          IntV_shapeI_gradshapeJ[i][j] +=
            fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeGrad(j, qp) * fe_vol_data.JxW(qp);

          IntV_shapeI_shapeJ[i][j] +=
            fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp) * fe_vol_data.JxW(qp);
        } // for qp

    cell_Gmatrices.push_back(std::move(IntV_shapeI_gradshapeJ));
    cell_Mmatrices.push_back(std::move(IntV_shapeI_shapeJ));

    const size_t num_faces = cell.faces_.size();
    VecMatDbl faces_Mmatrices;
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
      MatDbl IntS_shapeI_shapeJ(num_nodes, VecDbl(num_nodes, 0.0));
      for (unsigned int i = 0; i < num_nodes; ++i)
        for (unsigned int j = 0; j < num_nodes; ++j)
          for (const auto& qp : fe_srf_data.QuadraturePointIndices())
            IntS_shapeI_shapeJ[i][j] +=
              fe_srf_data.ShapeValue(i, qp) * fe_srf_data.ShapeValue(j, qp) * fe_srf_data.JxW(qp);

      faces_Mmatrices.push_back(std::move(IntS_shapeI_shapeJ));
    } // for face f

    cell_faceMmatrices.push_back(std::move(faces_Mmatrices));

  } // for cell

  opensn::log.Log() << "End cell matrices." << std::endl;

  // Make Grid internal face
  // mapping
  const auto cell_adj_mapping = sdm.MakeInternalFaceNodeMappings();

  // Define sweep chunk
  auto SweepChunk = [&ijk_mapping,
                     &grid,
                     &sdm,
                     &num_moments,
                     &phi_uk_man,
                     &psi_uk_man,
                     &m2d,
                     &d2m,
                     &phi_new,
                     &source_moments,
                     &psi_old,
                     &cell_Gmatrices,
                     &cell_Mmatrices,
                     &cell_faceMmatrices,
                     &cell_adj_mapping](const std::array<int64_t, 3>& ijk,
                                        const Vector3& omega,
                                        const size_t d,
                                        const MultiGroupXS& cell_xs)
  {
    const auto cell_global_id = ijk_mapping.MapNDtoLin(ijk[1], ijk[0], ijk[2]);
    const auto& cell = grid.cells[cell_global_id];
    const auto cell_local_id = cell.local_id_;
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const size_t num_faces = cell.faces_.size();

    const std::vector<double> zero_vector(num_groups, 0.0);

    const auto& G = cell_Gmatrices[cell_local_id];
    const auto& M = cell_Mmatrices[cell_local_id];

    MatDbl A(num_nodes, VecDbl(num_nodes, 0.0));
    MatDbl b(num_groups, VecDbl(num_nodes, 0.0));

    // Gradient matrix
    for (size_t i = 0; i < num_nodes; ++i)
      for (size_t j = 0; j < num_nodes; ++j)
        A[i][j] = omega.Dot(G[i][j]);

    // Surface integrals
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces_[f];
      const double mu = omega.Dot(face.normal_);

      if (mu < 0.0)
      {
        const auto& M_surf = cell_faceMmatrices[cell_local_id][f];

        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          for (size_t fj = 0; fj < num_face_nodes; ++fj)
          {
            const int j = cell_mapping.MapFaceNode(f, fj);

            const double* upwind_psi = zero_vector.data();
            if (face.has_neighbor_)
            {
              const auto& adj_cell = grid.cells[face.neighbor_id_];
              const int aj = cell_adj_mapping[cell.local_id_][f][fj];
              const int64_t ajmap = sdm.MapDOFLocal(adj_cell, aj, psi_uk_man, d, 0);
              upwind_psi = &psi_old[ajmap];
            }

            const double mu_Nij = -mu * M_surf[i][j];
            A[i][j] += mu_Nij;
            for (int g = 0; g < num_groups; ++g)
              b[g][i] += upwind_psi[g] * mu_Nij;
          } // for fj
        }   // for fi
      }     // if internal incident face
    }       // for face

    const auto& sigma_t = cell_xs.SigmaTotal();
    for (size_t g = 0; g < num_groups; ++g)
    {
      auto Atemp = A;
      VecDbl source(num_nodes, 0.0);
      // Nodal source moments
      for (size_t i = 0; i < num_nodes; ++i)
      {
        double temp_src = 0.0;
        for (size_t m = 0; m < num_moments; ++m)
        {
          const int64_t dof_map = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
          temp_src += m2d[m][d] * source_moments[dof_map];
        } // for m
        source[i] = temp_src;
      } // for i

      // Mass Matrix and Source
      const double sigma_tg = sigma_t[g];
      for (int i = 0; i < num_nodes; ++i)
      {
        double temp = 0.0;
        for (int j = 0; j < num_nodes; ++j)
        {
          const double Mij = M[i][j];
          Atemp[i][j] = A[i][j] + Mij * sigma_tg;
          temp += Mij * source[j];
        } // for j
        b[g][i] += temp;
      } // for i

      // Solve system
      GaussElimination(Atemp, b[g], static_cast<int>(num_nodes));
    } // for g

    // Accumulate flux-moments
    for (size_t m = 0; m < num_moments; ++m)
    {
      const double wn_d2m = d2m[m][d];
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell, i, phi_uk_man, m, 0);
        for (size_t g = 0; g < num_groups; ++g)
          phi_new[dof_map + g] += wn_d2m * b[g][i];
      }
    }

    // Save angular fluxes
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dof_map = sdm.MapDOFLocal(cell, i, psi_uk_man, d, 0);
      for (size_t g = 0; g < num_groups; ++g)
        psi_old[dof_map + g] = b[g][i];
    }
  };

  // Define sweep for all dirs
  auto Sweep = [&num_dirs, &quadrature, Nx, Ny, Nz, &SweepChunk, &xs]()
  {
    for (size_t d = 0; d < num_dirs; ++d)
    {
      const auto& omega = quadrature->omegas_[d];
      const auto& weight = quadrature->weights_[d];

      std::vector<int64_t> iorder, jorder, korder;
      if (omega.x > 0.0)
        iorder = Range<int64_t>(0, Nx);
      else
        iorder = Range<int64_t>(Nx - 1, -1, -1);
      if (omega.y > 0.0)
        jorder = Range<int64_t>(0, Ny);
      else
        jorder = Range<int64_t>(Ny - 1, -1, -1);
      if (omega.z > 0.0)
        korder = Range<int64_t>(0, Nz);
      else
        korder = Range<int64_t>(Nz - 1, -1, -1);

      for (auto i : iorder)
        for (auto j : jorder)
          for (auto k : korder)
            SweepChunk({i, j, k}, omega, d, xs);
    } // for d
  };

  // Define SetSource routine
  auto SetSource = [&source_moments,
                    &phi_old,
                    &q_source,
                    &grid,
                    &sdm,
                    &m_ell_em_map,
                    &xs,
                    num_moments,
                    &phi_uk_man]()
  {
    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();
      const auto& S = xs.TransferMatrices();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        for (size_t m = 0; m < num_moments; ++m)
        {
          const int64_t dof_map = sdm.MapDOFLocal(cell, i, phi_uk_man, m, 0);
          const auto ell = m_ell_em_map[m].ell;

          for (size_t g = 0; g < num_groups; ++g)
          {
            // Fixed source
            source_moments[dof_map + g] = q_source[dof_map + g];

            // Inscattering
            if (ell < S.size())
            {
              double inscat_g = 0.0;
              for (const auto& [row_g, gprime, sigma_sm] : S[ell].Row(g))
                inscat_g += sigma_sm * phi_old[dof_map + gprime];

              source_moments[dof_map + g] += inscat_g;
            }
          } // for g
        }   // for m
      }     // for node i
    }       // for cell
  };

  // Define L-infinite-norm
  auto ComputeRelativePWChange =
    [&grid, &sdm, &num_moments, &phi_uk_man](const std::vector<double>& in_phi_new,
                                             const std::vector<double>& in_phi_old)
  {
    double pw_change = 0.0;

    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        // Get scalar moments
        const int64_t m0_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

        const double* phi_new_m0 = &in_phi_new[m0_map];
        const double* phi_old_m0 = &in_phi_old[m0_map];
        for (size_t m = 0; m < num_moments; ++m)
        {
          const int64_t m_map = sdm.MapDOFLocal(cell, i, phi_uk_man, m, 0);

          const double* phi_new_m = &in_phi_new[m_map];
          const double* phi_old_m = &in_phi_old[m_map];

          for (size_t g = 0; g < num_groups; ++g)
          {
            const double abs_phi_new_g_m0 = std::fabs(phi_new_m0[g]);
            const double abs_phi_old_g_m0 = std::fabs(phi_old_m0[g]);

            const double max_denominator = std::max(abs_phi_new_g_m0, abs_phi_old_g_m0);

            const double delta_phi = std::fabs(phi_new_m[g] - phi_old_m[g]);

            if (max_denominator >= std::numeric_limits<double>::min())
              pw_change = std::max(delta_phi / max_denominator, pw_change);
            else
              pw_change = std::max(delta_phi, pw_change);
          } // for g
        }   // for m
      }     // for i
    }       // for cell

    return pw_change;
  };

  // Classic Richardson iteration
  opensn::log.Log() << "Starting iterations" << std::endl;
  for (size_t iter = 0; iter < 200; ++iter)
  {
    phi_new.assign(phi_new.size(), 0.0);
    // Build rhs
    SetSource();
    Sweep();

    const double rel_change = ComputeRelativePWChange(phi_new, phi_old);

    std::stringstream outstr;
    outstr << "Iteration " << std::setw(5) << iter << " ";
    {
      char buffer[100];
      snprintf(buffer, 100, "%11.3e\n", rel_change);
      outstr << buffer;
    }

    opensn::log.Log() << outstr.str();

    phi_old = phi_new;

    if (rel_change < 1.0e-6 and iter > 0)
      break;
  } // for iteration

  // Create Field Functions
  std::vector<std::shared_ptr<FieldFunctionGridBased>> ff_list;

  ff_list.push_back(std::make_shared<FieldFunctionGridBased>(
    "Phi", sdm_ptr, Unknown(UnknownType::VECTOR_N, num_groups)));

  const std::vector<std::string> dim_strings = {"x", "y", "z"};
  for (const std::string& dim : dim_strings)
    ff_list.push_back(std::make_shared<FieldFunctionGridBased>(
      "J-" + dim, sdm_ptr, Unknown(UnknownType::VECTOR_N, num_groups)));

  // Localize zeroth moment
  // This routine extracts a single moment vector
  // from the vector that contains multiple moments
  const UnknownManager m0_uk_man({Unknown(UnknownType::VECTOR_N, num_groups)});
  const size_t num_m0_dofs = sdm.GetNumLocalDOFs(m0_uk_man);

  std::vector<double> m0_phi(num_m0_dofs, 0.0);
  std::vector<double> mx_phi(num_m0_dofs, 0.0); // Y(1,1)  - X-component
  std::vector<double> my_phi(num_m0_dofs, 0.0); // Y(1,-1) - Y-component
  std::vector<double> mz_phi(num_m0_dofs, 0.0); // Y(1,0)  - Z-component

  sdm.CopyVectorWithUnknownScope(phi_old, m0_phi, phi_uk_man, 0, m0_uk_man, 0);

  ff_list[0]->UpdateFieldVector(m0_phi);

  std::array<unsigned int, 3> j_map = {0, 0, 0};
  if (dimension == 1 and num_moments >= 2)
    j_map = {0, 0, 1};
  if (dimension == 2 and num_moments >= 3)
    j_map = {2, 1, 0};
  if (dimension == 3 and num_moments >= 4)
    j_map = {3, 1, 2};

  sdm.CopyVectorWithUnknownScope(phi_old, mx_phi, phi_uk_man, j_map[0], m0_uk_man, 0);
  sdm.CopyVectorWithUnknownScope(phi_old, my_phi, phi_uk_man, j_map[1], m0_uk_man, 0);
  sdm.CopyVectorWithUnknownScope(phi_old, mz_phi, phi_uk_man, j_map[2], m0_uk_man, 0);

  ff_list[1]->UpdateFieldVector(mx_phi);
  ff_list[2]->UpdateFieldVector(my_phi);
  ff_list[3]->UpdateFieldVector(mz_phi);

  // Update field function
  FieldFunctionGridBased::FFList const_ff_list;
  for (const auto& ff_ptr : ff_list)
    const_ff_list.push_back(ff_ptr);
  FieldFunctionGridBased::ExportMultipleToVTK("SimTest_91a_PWLD", const_ff_list);

  return ParameterBlock();
}

} // namespace unit_sim_tests
