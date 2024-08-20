#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/random_number_generation/random_number_generator.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "lua/framework/console/console.h"
#include "framework/mesh/raytrace/raytracer.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

using namespace opensn;

namespace unit_sim_tests
{

ParameterBlock SimTest93_RayTracing(const InputParameters& params);

RegisterWrapperFunctionInNamespace(unit_tests, SimTest93_RayTracing, nullptr, SimTest93_RayTracing);

ParameterBlock
SimTest93_RayTracing(const InputParameters&)
{
  const std::string fname = "SimTest93_RayTracing";
  opensn::log.Log() << "SimTest93_RayTracing";

  // Get grid
  auto grid_ptr = GetCurrentMesh();
  const auto& grid = *grid_ptr;

  opensn::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

  const auto dimension = grid.Dimension();

  // Set parameters
  const size_t num_groups = 1;
  const size_t scattering_order = 1;
  const auto& L = scattering_order;
  const size_t num_moments = (dimension == 1)   ? L + 1
                             : (dimension == 2) ? (L + 1) * (L + 2) / 2
                             : (dimension == 3) ? (L + 1) * (L + 1)
                                                : 0;
  const double sigma_t = 0.27;

  // Build harmonic map
  std::vector<std::pair<int, int>> m_to_ell_em_map;
  if (dimension == 1)
    for (int ell = 0; ell <= scattering_order; ++ell)
      m_to_ell_em_map.emplace_back(ell, 0);
  else if (dimension == 2)
    for (int ell = 0; ell <= scattering_order; ++ell)
      for (int m = -ell; m <= ell; m += 2)
        m_to_ell_em_map.emplace_back(ell, m);
  else if (dimension == 3)
    for (int ell = 0; ell <= scattering_order; ++ell)
      for (int m = -ell; m <= ell; ++m)
        m_to_ell_em_map.emplace_back(ell, m);

  // Make SDM
  std::shared_ptr<SpatialDiscretization> sdm_ptr = PieceWiseLinearDiscontinuous::New(grid);
  const auto& sdm = *sdm_ptr;

  UnknownManager phi_uk_man;
  for (size_t m = 0; m < num_moments; ++m)
    phi_uk_man.AddUnknown(UnknownType::VECTOR_N, num_groups);

  const size_t num_fem_local_dofs = sdm.GetNumLocalDOFs(phi_uk_man);
  const size_t num_fem_globl_dofs = sdm.GetNumGlobalDOFs(phi_uk_man);

  opensn::log.Log() << "Num local FEM DOFs: " << num_fem_local_dofs;
  opensn::log.Log() << "Num globl FEM DOFs: " << num_fem_globl_dofs;

  // Define tallies
  std::vector<double> phi_tally(num_fem_local_dofs, 0.0);

  // Define particle data structure
  struct Particle
  {
    Vector3 position = {0.0, 0.0, 0.0};
    Vector3 direction = {0.0, 0.0, 0.0};
    int energy_group = 0;
    double weight = 1.0;

    uint64_t cell_id = 0;

    bool alive = true;
  };

  // Define source position
  //                                              and find cell containing it
  const Vector3 source_pos = {0.0, 0.0, 0.0};

  Cell const* source_cell_ptr = nullptr;

  for (auto& cell : grid.local_cells)
    if (grid.CheckPointInsideCell(cell, source_pos))
    {
      source_cell_ptr = &cell;
      break;
    }
  if (source_cell_ptr == nullptr)
    throw std::logic_error(fname + ": Source cell not found.");

  const uint64_t source_cell_id = source_cell_ptr->global_id;

  // Define lambdas
  RandomNumberGenerator rng;
  auto SampleRandomDirection = [&rng]()
  {
    double costheta = 2.0 * rng.Rand() - 1.0;
    double theta = acos(costheta);
    double varphi = rng.Rand() * 2.0 * M_PI;

    return Vector3{sin(theta) * cos(varphi), sin(theta) * sin(varphi), cos(theta)};
  };

  auto ContributePWLDTally =
    [&sdm, &grid, &phi_tally, &phi_uk_man, &sigma_t, &num_moments, &m_to_ell_em_map](
      const Cell& cell,
      const Vector3& positionA,
      const Vector3& positionB,
      const Vector3& omega,
      const int g,
      double weight)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto phi_theta = OmegaToPhiThetaSafe(omega);
    const double phi = phi_theta.first;
    const double theta = phi_theta.second;

    std::vector<double> segment_lengths;
    PopulateRaySegmentLengths(grid, cell, positionA, positionB, omega, segment_lengths);

    Vector<double> shape_values_k;   // At s_k
    Vector<double> shape_values_kp1; // At s_{k+1}

    cell_mapping.ShapeValues(positionA, shape_values_k);

    double d_run_sum = 0.0;
    for (const auto& segment_length_k : segment_lengths)
    {
      d_run_sum += segment_length_k;
      const double& d = d_run_sum;

      cell_mapping.ShapeValues(positionA + omega * d, shape_values_kp1);

      const auto& b_ik = shape_values_k;
      const auto& b_ikp1 = shape_values_kp1;
      const double& ell_k = segment_length_k;

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const double C0 = b_ik(i) * ell_k;
        const double C1 = b_ikp1(i) - b_ik(i);

        for (size_t m = 0; m < num_moments; ++m)
        {
          const int64_t dof_map = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);

          // Apply harmonic weight
          const auto& ell_em = m_to_ell_em_map.at(m);
          const int ell = ell_em.first;
          const int em = ell_em.second;

          double w_harmonic = Ylm(ell, em, phi, theta);

          // Apply exponential attenuation weight
          double w_exp =
            (C0 / sigma_t) * (1.0 - exp(-sigma_t * ell_k)) +
            (C1 / (sigma_t * sigma_t)) * (1.0 - (1 + sigma_t * ell_k) * exp(-sigma_t * ell_k));
          w_exp *= weight / (ell_k * ell_k);

          // Combine
          double w_avg = w_harmonic * w_exp;

          phi_tally[dof_map] += ell_k * w_avg;
        } // for moment m
      }   // for node i

      shape_values_k = shape_values_kp1;
      weight *= exp(-sigma_t * segment_length_k);
    } // for d
  };

  auto GetCellApproximateSize = [&grid](const Cell& cell)
  {
    const auto& v0 = grid.vertices[cell.vertex_ids[0]];
    double xmin = v0.x, xmax = v0.x;
    double ymin = v0.y, ymax = v0.y;
    double zmin = v0.z, zmax = v0.z;

    for (uint64_t vid : cell.vertex_ids)
    {
      const auto& v = grid.vertices[vid];

      xmin = std::min(xmin, v.x);
      xmax = std::max(xmax, v.x);
      ymin = std::min(ymin, v.y);
      ymax = std::max(ymax, v.y);
      zmin = std::min(zmin, v.z);
      zmax = std::max(zmax, v.z);
    }

    return (Vector3(xmin, ymin, zmin) - Vector3(xmax, ymax, zmax)).Norm();
  };

  // Create raytracer
  std::vector<double> cell_sizes(grid.local_cells.size(), 0.0);
  for (const auto& cell : grid.local_cells)
    cell_sizes[cell.local_id] = GetCellApproximateSize(cell);

  RayTracer ray_tracer(grid, cell_sizes);

  // Run rays
  const auto PWLD = SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;

  const size_t num_particles = 100'000;
  for (size_t n = 0; n < num_particles; ++n)
  {
    if (n % size_t(num_particles / 10.0) == 0)
      std::cout << "#particles = " << n << "\n";
    // Create the particle
    const auto omega = SampleRandomDirection();
    Particle particle{source_pos, omega, 0, 1.0, source_cell_id, true};

    while (particle.alive)
    {
      // Get the current cell
      const auto& cell = grid.cells[particle.cell_id];

      // Perform the trace to the next surface
      auto destination_info = ray_tracer.TraceRay(cell, particle.position, particle.direction);

      const Vector3& end_of_track_position = destination_info.pos_f;

      // Make tally contribution
      ContributePWLDTally(cell,
                          particle.position,
                          end_of_track_position,
                          particle.direction,
                          particle.energy_group,
                          particle.weight);

      // Process cell transfer
      //                                or death
      if (not destination_info.particle_lost)
      {
        const auto& f = destination_info.destination_face_index;
        const auto& current_cell_face = cell.faces[f];

        if (current_cell_face.has_neighbor)
          particle.cell_id = current_cell_face.neighbor_id;
        else
          particle.alive = false; // Death at the boundary
      }
      else
      {
        std::cout << "particle" << n << " lost " << particle.position.PrintStr() << " "
                  << particle.direction.PrintStr() << " "
                  << "\n";
        break;
      }

      const auto& pA = particle.position;
      const auto& pB = end_of_track_position;
      particle.weight *= exp(-sigma_t * (pB - pA).Norm()); // Attenuation
      particle.position = end_of_track_position;
    } // while ray alive

  } // for ray n

  // Post process tallies
  for (const auto& cell : grid.local_cells)
  {
    // Compute mass matrix and its inverse
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
    const size_t num_nodes = cell_mapping.NumNodes();

    DenseMatrix<double> M(num_nodes, num_nodes, 0.0);
    for (auto qp : fe_vol_data.QuadraturePointIndices())
      for (size_t i = 0; i < num_nodes; ++i)
        for (size_t j = 0; j < num_nodes; ++j)
          M(i, j) +=
            fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp) * fe_vol_data.JxW(qp);

    auto M_inv = Inverse(M);

    // Apply projection
    Vector<double> T(num_nodes, 0.0);
    for (size_t m = 0; m < num_moments; ++m)
      for (size_t g = 0; g < num_groups; ++g)
      {
        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
          T(i) = phi_tally[imap] / num_particles;
        }

        auto phi_uc = MatMul(M_inv, T);

        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imap = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
          phi_tally[imap] = phi_uc(i);
        }
      } // for group g

  } // for cell

  // Create Field Functions
  std::vector<std::shared_ptr<FieldFunctionGridBased>> ff_list;

  ff_list.push_back(std::make_shared<FieldFunctionGridBased>(
    "Phi", sdm_ptr, Unknown(UnknownType::VECTOR_N, num_groups)));

  // Localize zeroth moment
  // This routine extracts a single moment vector
  // from the vector that contains multiple moments
  const UnknownManager m0_uk_man({Unknown(UnknownType::VECTOR_N, num_groups)});
  const size_t num_m0_dofs = sdm.GetNumLocalDOFs(m0_uk_man);

  std::vector<double> m0_phi(num_m0_dofs, 0.0);

  sdm.CopyVectorWithUnknownScope(phi_tally, m0_phi, phi_uk_man, 0, m0_uk_man, 0);

  ff_list[0]->UpdateFieldVector(m0_phi);

  // Update field function
  std::vector<std::shared_ptr<const FieldFunctionGridBased>> const_ff_list;
  for (const auto& ff_ptr : ff_list)
    const_ff_list.push_back(ff_ptr);
  FieldFunctionGridBased::ExportMultipleToVTK("SimTest_93", const_ff_list);

  return ParameterBlock();
}

} // namespace unit_sim_tests
