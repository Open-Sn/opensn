#include "framework/math/spatial_discretization/finite_volume/finite_volume.h"
#include "framework/physics/physics_material/multi_group_xs/multi_group_xs.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/math_range.h"
#include "framework/data_types/ndarray.h"
#include "lua/framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>

using namespace opensn;

namespace unit_sim_tests
{
typedef AngularQuadrature::HarmonicIndices YlmIndices;
double ComputeRelativePWChange(const MeshContinuum& grid,
                               const SpatialDiscretization& sdm,
                               const UnknownManager& phi_uk_man,
                               const std::vector<double>& in_phi_new,
                               const std::vector<double>& in_phi_old);

std::vector<double> SetSource(const MeshContinuum& grid,
                              const SpatialDiscretization& sdm,
                              const UnknownManager& phi_uk_man,
                              const std::vector<double>& q_source,
                              const std::vector<double>& phi_old,
                              const MGXS& xs,
                              const std::vector<YlmIndices>& m_ell_em_map);

/**WDD Sweep. */
ParameterBlock SimTest06_WDD(const InputParameters&);

RegisterWrapperFunctionNamespace(unit_testsB, SimTest06_WDD, nullptr, SimTest06_WDD);

ParameterBlock
SimTest06_WDD(const InputParameters&)
{
  const std::string fname = "SimTest06_WDD";

  opensn::log.Log() << "SimTest06_WDD num_args = " << 0;

  if (opensn::mpi_comm.size() != 1)
    throw std::logic_error(fname + ": Is serial only.");

  // Get grid
  auto grid_ptr = GetCurrentMesh();
  const auto& grid = *grid_ptr;

  opensn::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  // Make Orthogonal mapping
  const auto ijk_info = grid.GetIJKInfo();
  const auto& ijk_mapping = grid.MakeIJKToGlobalIDMapping();
  const auto cell_ortho_sizes = grid.MakeCellOrthoSizes();

  const auto Nx = static_cast<int64_t>(ijk_info[0]);
  const auto Ny = static_cast<int64_t>(ijk_info[1]);
  const auto Nz = static_cast<int64_t>(ijk_info[2]);

  const auto Dim1 = DIMENSION_1;
  const auto Dim2 = DIMENSION_2;
  const auto Dim3 = DIMENSION_3;

  int dimension = 0;
  if (grid.Attributes() & Dim1)
    dimension = 1;
  if (grid.Attributes() & Dim2)
    dimension = 2;
  if (grid.Attributes() & Dim3)
    dimension = 3;

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = FiniteVolume::New(grid);
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
  MGXS xs;
  xs.Initialize("xs_graphite_pure.xs");

  // Initializes vectors
  std::vector<double> phi_old(num_local_phi_dofs, 0.0);
  std::vector<double> psi(num_local_psi_dofs, 0.0);
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
      } // for node i
    }   // if inside box
  }     // for cell

  // Define sweep chunk
  typedef NDArray<double> IJKArrayDbl;
  IJKArrayDbl psi_ds_x(std::array<int64_t, 4>{Nx, Ny, Nz, num_groups});
  IJKArrayDbl psi_ds_y(std::array<int64_t, 4>{Nx, Ny, Nz, num_groups});
  IJKArrayDbl psi_ds_z(std::array<int64_t, 4>{Nx, Ny, Nz, num_groups});

  typedef Vector3 Vec3;
  auto SweepChunk = [&ijk_info,
                     &ijk_mapping,
                     &cell_ortho_sizes,
                     &grid,
                     &sdm,
                     &num_moments,
                     &phi_uk_man,
                     &psi_uk_man,
                     &m2d,
                     &d2m,
                     &phi_new,
                     &source_moments,
                     &psi,
                     &psi_ds_x,
                     &psi_ds_y,
                     &psi_ds_z](const std::array<int64_t, 3>& ijk,
                                const Vec3& omega,
                                const size_t d,
                                const MGXS& cell_xs)
  {
    const auto cell_global_id = ijk_mapping.MapNDtoLin(ijk[1], ijk[0], ijk[2]);
    const auto& cell = grid.cells[cell_global_id];
    const auto cell_local_id = cell.local_id_;

    const auto& cell_ortho_size = cell_ortho_sizes[cell_local_id];
    const double dx = cell_ortho_size.x;
    const double dy = cell_ortho_size.y;
    const double dz = cell_ortho_size.z;

    const auto i = ijk[0];
    const auto Nx = ijk_info[0];
    const auto j = ijk[1];
    const auto Ny = ijk_info[1];
    const auto k = ijk[2];
    const auto Nz = ijk_info[2];

    const std::vector<double> zero_vector(num_groups, 0.0);

    const double* psi_us_x = zero_vector.data();
    const double* psi_us_y = zero_vector.data();
    const double* psi_us_z = zero_vector.data();

    if (omega.x > 0.0 and i > 0)
      psi_us_x = &psi_ds_x(i - 1, j, k, 0);
    if (omega.x < 0.0 and i < (Nx - 1))
      psi_us_x = &psi_ds_x(i + 1, j, k, 0);
    if (omega.y > 0.0 and j > 0)
      psi_us_y = &psi_ds_y(i, j - 1, k, 0);
    if (omega.y < 0.0 and j < (Ny - 1))
      psi_us_y = &psi_ds_y(i, j + 1, k, 0);
    if (omega.z > 0.0 and k > 0)
      psi_us_z = &psi_ds_z(i, j, k - 1, 0);
    if (omega.z < 0.0 and k < (Nz - 1))
      psi_us_z = &psi_ds_z(i, j, k + 1, 0);

    const auto& sigma_t = cell_xs.SigmaTotal();
    for (size_t g = 0; g < num_groups; ++g)
    {
      double rhs = 0.0;
      // Source moments
      for (size_t m = 0; m < num_moments; ++m)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell, 0, phi_uk_man, m, g);
        rhs += source_moments[dof_map] * m2d[m][d];
      }

      if (Nx > 1)
        rhs += 2.0 * std::fabs(omega.x) * psi_us_x[g] / dx;
      if (Ny > 1)
        rhs += 2.0 * std::fabs(omega.y) * psi_us_y[g] / dy;
      if (Nz > 1)
        rhs += 2.0 * std::fabs(omega.z) * psi_us_z[g] / dz;

      double lhs = sigma_t[g];
      if (Nx > 1)
        lhs += 2.0 * std::fabs(omega.x) / dx;
      if (Ny > 1)
        lhs += 2.0 * std::fabs(omega.y) / dy;
      if (Nz > 1)
        lhs += 2.0 * std::fabs(omega.z) / dz;

      double psi_ijk = rhs / lhs;

      // Accumulate flux-moments
      for (size_t m = 0; m < num_moments; ++m)
      {
        const int64_t dof_map = sdm.MapDOFLocal(cell, 0, phi_uk_man, m, g);
        phi_new[dof_map] += d2m[m][d] * psi_ijk;
      }

      // Save angular fluxes
      const int64_t psi_map = sdm.MapDOFLocal(cell, 0, psi_uk_man, d, g);
      psi[psi_map] = psi_ijk;

      psi_ds_x(i, j, k, g) = 2.0 * psi_ijk - psi_us_x[g];
      psi_ds_y(i, j, k, g) = 2.0 * psi_ijk - psi_us_y[g];
      psi_ds_z(i, j, k, g) = 2.0 * psi_ijk - psi_us_z[g];
    } // for g
  };

  // Define sweep for all dirs
  auto Sweep = [&num_dirs, &quadrature, Nx, Ny, Nz, &SweepChunk, &xs]()
  {
    for (size_t d = 0; d < num_dirs; ++d)
    {
      const auto& omega = quadrature->omegas_[d];

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

  // Classic Richardson iteration
  opensn::log.Log() << "Starting iterations" << std::endl;
  for (size_t iter = 0; iter < 200; ++iter)
  {
    phi_new.assign(phi_new.size(), 0.0);
    // Build rhs
    source_moments = SetSource(grid, sdm, phi_uk_man, q_source, phi_old, xs, m_ell_em_map);
    Sweep();

    const double rel_change = ComputeRelativePWChange(grid, sdm, phi_uk_man, phi_new, phi_old);

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

  // Localize zeroth moment
  // This routine extracts a single moment vector
  // from the vector that contains multiple moments
  const UnknownManager m0_uk_man({Unknown(UnknownType::VECTOR_N, num_groups)});
  const size_t num_m0_dofs = sdm.GetNumLocalDOFs(m0_uk_man);

  std::vector<double> m0_phi(num_m0_dofs, 0.0);

  sdm.CopyVectorWithUnknownScope(phi_old, m0_phi, phi_uk_man, 0, m0_uk_man, 0);

  // Create Field Function
  auto phi_ff = std::make_shared<FieldFunctionGridBased>(
    "Phi", sdm_ptr, Unknown(UnknownType::VECTOR_N, num_groups));

  phi_ff->UpdateFieldVector(m0_phi);

  FieldFunctionGridBased::ExportMultipleToVTK("SimTest_06_WDD", {phi_ff});

  return ParameterBlock();
}

double
ComputeRelativePWChange(const MeshContinuum& grid,
                        const SpatialDiscretization& sdm,
                        const UnknownManager& phi_uk_man,
                        const std::vector<double>& in_phi_new,
                        const std::vector<double>& in_phi_old)
{
  double pw_change = 0.0;
  const size_t num_moments = phi_uk_man.unknowns_.size();
  const size_t num_groups = phi_uk_man.unknowns_.front().num_components_;

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
}

std::vector<double>
SetSource(const MeshContinuum& grid,
          const SpatialDiscretization& sdm,
          const UnknownManager& phi_uk_man,
          const std::vector<double>& q_source,
          const std::vector<double>& phi_old,
          const MGXS& xs,
          const std::vector<YlmIndices>& m_ell_em_map)
{
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(phi_uk_man);
  std::vector<double> source_moments(num_local_phi_dofs, 0.0);

  const size_t num_moments = phi_uk_man.unknowns_.size();
  const size_t num_groups = phi_uk_man.unknowns_.front().num_components_;

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

  return source_moments;
}

} // namespace unit_sim_tests
