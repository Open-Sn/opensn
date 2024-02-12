#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/lbs_adj_solver.h"
#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/lbs_adjoint.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/source_functions/source_function.h"

#include "framework/physics/physics_material/multi_group_xs/adjoint_mgxs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/math_vector_nx.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include <utility>
#include <fstream>

namespace opensn::lbs
{

OpenSnRegisterObjectInNamespace(lbs, DiscreteOrdinatesAdjointSolver);

InputParameters
DiscreteOrdinatesAdjointSolver::GetInputParameters()
{
  InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();

  params.SetGeneralDescription("Adjoint capability");

  params.SetClassName("DiscreteOrdinatesAdjointSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesAdjointSolver");

  return params;
}

DiscreteOrdinatesAdjointSolver::DiscreteOrdinatesAdjointSolver(const InputParameters& params)
  : lbs::DiscreteOrdinatesSolver(params)
{
}

DiscreteOrdinatesAdjointSolver::DiscreteOrdinatesAdjointSolver(const std::string& solver_name)
  : lbs::DiscreteOrdinatesSolver(solver_name)
{
}

void
DiscreteOrdinatesAdjointSolver::Initialize()
{
  LBSSolver::Initialize();

  // Create adjoint cross sections
  std::map<int, std::shared_ptr<MultiGroupXS>> matid_to_adj_xs_map;
  for (const auto& [matid, xs] : matid_to_xs_map_)
    matid_to_adj_xs_map[matid] =
      std::make_shared<AdjointMGXS>(*std::dynamic_pointer_cast<MultiGroupXS>(xs));
  matid_to_xs_map_ = std::move(matid_to_adj_xs_map);

  for (const auto& cell : grid_ptr_->local_cells)
    cell_transport_views_[cell.local_id_].ReassignXS(*matid_to_xs_map_[cell.material_id_]);

  // Initialize source func
  using namespace std::placeholders;
  auto src_function = std::make_shared<SourceFunction>(*this);
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);

  // Initialize groupsets for sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }
  InitializeSolverSchemes();

  source_event_tag_ = log.GetRepeatingEventTag("Set Source");
}

void
DiscreteOrdinatesAdjointSolver::Execute()
{
  primary_ags_solver_->Setup();
  primary_ags_solver_->Solve();

  // Apply post-processing
  log.Log() << "LBAdjointSolver: post-processing.";
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets_)
    for (const auto& group : groupset.groups_)
      set_group_numbers.insert(group.id_);

  const auto& m_to_ell_em_map = groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap();

  // Reorient phi-moments for reverse angle
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_view = cell_transport_views_[cell.local_id_];
    const int num_nodes = cell_view.NumNodes();

    for (int i = 0; i < num_nodes; ++i)
    {
      for (int m = 0; m < num_moments_; ++m)
      {
        const auto& ell = m_to_ell_em_map[m].ell;
        const auto dof_map = cell_view.MapDOF(i, m, 0);

        for (int g : set_group_numbers)
          phi_old_local_[dof_map + g] *= pow(-1.0, ell);
      } // for moment
    }   // node i
  }     // for cell

  UpdateFieldFunctions();
}

double
DiscreteOrdinatesAdjointSolver::ComputeInnerProduct()
{
  double local_integral = 0.0;

  // Material sources
  for (const auto& cell : grid_ptr_->local_cells)
  {
    if (matid_to_src_map_.count(cell.material_id_) == 0)
      continue; // Skip if no src

    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& source = matid_to_src_map_[cell.material_id_];
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];

    for (const auto& group : groups_)
    {
      const auto& g = group.id_;
      const auto& q = source->source_value_g_[g];

      if (q > 0.0)
      {
        for (int i = 0; i < transport_view.NumNodes(); ++i)
        {
          const auto dof_map = transport_view.MapDOF(i, 0, g);
          const auto& phi = phi_old_local_[dof_map];

          local_integral += q * phi * fe_values.intV_shapeI[i];
        } // for node
      }   // check source value >0
    }     // for group
  }       // for cell

  // Point sources
  for (const auto& point_source : point_sources_)
  {
    for (const auto& subscriber : point_source.Subscribers())
    {
      const auto& cell = grid_ptr_->local_cells[subscriber.cell_local_id];
      const auto& transport_view = cell_transport_views_[cell.local_id_];
      const auto& source_strength = point_source.Strength();
      const auto& shape_values = subscriber.shape_values;

      for (const auto& group : groups_)
      {
        const auto& g = group.id_;
        const auto& S = source_strength[g] * subscriber.volume_weight;

        if (S > 0.0)
        {
          const auto num_nodes = transport_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            const auto dof_map = transport_view.MapDOF(i, 0, g);
            const auto& phi = phi_old_local_[dof_map];

            local_integral += S * phi * shape_values[i];
          } // for node
        }   // check source value >0
      }     // for group
    }       // for cell
  }         // for point source

  // Distributed sources
  for (const auto& distributed_source : distributed_sources_)
  {
    for (const auto& local_id : distributed_source.Subscribers())
    {
      const auto& cell = grid_ptr_->local_cells[local_id];
      const auto& transport_view = cell_transport_views_[local_id];
      const auto& fe_values = unit_cell_matrices_[local_id];
      const auto nodes = discretization_->GetCellNodeLocations(cell);

      for (int i = 0; i < transport_view.NumNodes(); ++i)
      {
        // Compute group-wise values for this node
        const auto src = distributed_source(cell, nodes[i], num_groups_);

        // Contribute to the source moments
        const auto& intV_shapeI = fe_values.intV_shapeI[i];
        const auto dof_map = transport_view.MapDOF(i, 0, 0);
        for (const auto& group : groups_)
          local_integral += src[group.id_] * intV_shapeI;
      }
    }
  }

  double global_integral = 0.0;
  mpi_comm.all_reduce(local_integral, global_integral, mpi::op::sum<double>());
  return global_integral;
}

void
DiscreteOrdinatesAdjointSolver::ExportImportanceMap(const std::string& file_name)
{
  // Determine cell averaged importance map
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets_)
    for (const auto& group : groupset.groups_)
      set_group_numbers.insert(group.id_);

  const auto& m_to_ell_em_map = groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap();

  typedef VectorN<4> Arr4; // phi, J_x, J_y, J_z
  typedef std::vector<Arr4> MGVec4;
  typedef std::vector<MGVec4> VecOfMGVec4;
  const size_t num_groups = set_group_numbers.size();
  const size_t num_cells = grid_ptr_->local_cells.size();

  VecOfMGVec4 cell_avg_p1_moments(num_cells, MGVec4(num_groups));
  {

    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_view = cell_transport_views_[cell.local_id_];
      const int num_nodes = cell_view.NumNodes();
      const auto& fe_values = unit_cell_matrices_[cell.local_id_];

      VecOfMGVec4 nodal_p1_moments(num_nodes);
      for (int i = 0; i < num_nodes; ++i)
      {
        // Get multigroup p1_moments
        MGVec4 p1_moments(set_group_numbers.size(), VecDbl{0.0, 0.0, 0.0, 0.0});
        for (int m = 0; m < std::max(static_cast<int>(num_moments_), 4); ++m)
        {
          const auto& ell = m_to_ell_em_map[m].ell;
          const auto& em = m_to_ell_em_map[m].m;

          size_t dof_map = cell_view.MapDOF(i, m, 0); // unknown map

          for (int g : set_group_numbers)
          {
            if (ell == 0 and em == 0)
              p1_moments[g](0) = std::fabs(phi_old_local_[dof_map + g]);
            if (ell == 1 and em == 1)
              p1_moments[g](1) = phi_old_local_[dof_map + g];
            if (ell == 1 and em == -1)
              p1_moments[g](2) = phi_old_local_[dof_map + g];
            if (ell == 1 and em == 0)
              p1_moments[g](3) = phi_old_local_[dof_map + g];
          } // for g
        }   // for m

        nodal_p1_moments[i] = std::move(p1_moments);
      } // for node i

      // Determine nodal average p1_moments
      for (int g : set_group_numbers)
      {
        VectorN<4> cell_p1_avg(VecDbl{0.0, 0.0, 0.0, 0.0});

        double volume_total = 0.0;
        for (int i = 0; i < num_nodes; ++i)
        {
          double IntV_shapeI = fe_values.intV_shapeI[i];
          cell_p1_avg += nodal_p1_moments[i][g] * IntV_shapeI;
          volume_total += IntV_shapeI;
        } // for node i
        cell_p1_avg /= volume_total;

        cell_avg_p1_moments[cell.local_id_][g] = cell_p1_avg;
      } // for g
    }   // for cell
  }

  // Determine cell-based exponential-representations
  typedef std::pair<double, double> ABCoeffsPair;
  typedef std::vector<ABCoeffsPair> VecOfABCoeffsPair;
  typedef std::vector<VecOfABCoeffsPair> ExpReps;
  ExpReps cell_exp_reps(num_cells, VecOfABCoeffsPair(num_groups, {0.0, 0.0}));
  {
    for (const auto& cell : grid_ptr_->local_cells)
    {
      for (int g : set_group_numbers)
      {
        const auto& p1_moments = cell_avg_p1_moments[cell.local_id_][g];

        auto a_b = MakeExpRepFromP1({p1_moments[0], p1_moments[1], p1_moments[2], p1_moments[3]});

        cell_exp_reps[cell.local_id_][g] = std::make_pair(a_b[0], a_b[1]);
      } // for g
    }   // for cell
  }

  log.Log() << "Exporting importance map to binary file " << file_name;

  const auto locJ_io_flags = std::ofstream::binary | std::ofstream::out;
  const auto loc0_io_flags = locJ_io_flags | std::ofstream::trunc;
  const bool is_home = (opensn::mpi_comm.rank() == 0);

  // Build header
  std::string header_info = "Chi-Tech LinearBoltzmann: Importance map file\n"
                            "Header size: 500 bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t num_global_cells\n"
                            "uint64_t num_groups\n"
                            "uint64_t num_records\n"
                            "Each record:\n"
                            "  uint64_t     cell_global_id\n"
                            "  unsigned int group_num\n"
                            "  double       phi\n"
                            "  double       J_x\n"
                            "  double       J_y\n"
                            "  double       J_z\n"
                            "  double       a_coefficient\n"
                            "  double       b_coefficient\n";

  int header_size = (int)header_info.length();

  char header_bytes[400];
  memset(header_bytes, '-', 400);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, 399));
  header_bytes[399] = '\0';

  // Process each location
  uint64_t num_global_cells = grid_ptr_->GetGlobalNumberOfCells();
  for (int locationJ = 0; locationJ < opensn::mpi_comm.size(); ++locationJ)
  {
    log.LogAll() << "  Barrier at " << locationJ;
    opensn::mpi_comm.barrier();
    if (opensn::mpi_comm.rank() != locationJ)
      continue;

    log.LogAll() << "  Location " << locationJ << " appending data.";

    std::ofstream file(file_name, is_home ? loc0_io_flags : locJ_io_flags);
    ChiLogicalErrorIf(not file.is_open(),
                      std::string(__FUNCTION__) + ": Location " + std::to_string(mpi_comm.rank()) +
                        ", failed to open file " + file_name + ".");

    if (is_home)
    {
      file << header_bytes;
      uint64_t num_groups_t = groups_.size();
      uint64_t num_records = num_global_cells * num_groups;

      file.write((char*)&num_global_cells, sizeof(uint64_t));
      file.write((char*)&num_groups_t, sizeof(uint64_t));
      file.write((char*)&num_records, sizeof(uint64_t));
    }

    for (const auto& cell : grid_ptr_->local_cells)
    {
      MGVec4 p1_moments = cell_avg_p1_moments[cell.local_id_];
      VecOfABCoeffsPair exp_rep = cell_exp_reps[cell.local_id_];

      auto cell_global_id = static_cast<uint64_t>(cell.global_id_);
      for (int group : set_group_numbers)
      {
        auto g = static_cast<unsigned int>(group);
        file.write((char*)&cell_global_id, sizeof(uint64_t));
        file.write((char*)&g, sizeof(unsigned int));

        for (int m = 0; m < 4; ++m)
          file.write((char*)&p1_moments[group](m), sizeof(double));

        file.write((char*)&exp_rep[group].first, sizeof(double));
        file.write((char*)&exp_rep[group].second, sizeof(double));
      } // for g
    }   // for cell

    file.close();
  } // for location

  log.LogAll() << "Done exporting importance map to binary file " << file_name;
  opensn::mpi_comm.barrier();
}

} // namespace opensn::lbs
