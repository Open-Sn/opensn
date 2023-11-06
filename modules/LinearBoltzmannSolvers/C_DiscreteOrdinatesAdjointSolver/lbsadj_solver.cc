#include "modules/LinearBoltzmannSolvers/C_DiscreteOrdinatesAdjointSolver/lbsadj_solver.h"
#include "modules/LinearBoltzmannSolvers/C_DiscreteOrdinatesAdjointSolver/lbs_adjoint.h"
#include "framework/ChiObjectFactory.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/SourceFunctions/adjoint_src_function.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "framework/physics/PhysicsMaterial/MultiGroupXS/adjoint_mgxs.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/mesh/LogicalVolume/LogicalVolume.h"
#include "framework/math/chi_math_vectorNX.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/mpi/chi_mpi.h"
#include <utility>
#include <fstream>

namespace lbs
{

RegisterChiObject(lbs, DiscreteOrdinatesAdjointSolver);

chi::InputParameters
DiscreteOrdinatesAdjointSolver::GetInputParameters()
{
  chi::InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();

  params.SetGeneralDescription("Adjoint capability");

  params.SetClassName("DiscreteOrdinatesAdjointSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesAdjointSolver");

  return params;
}

DiscreteOrdinatesAdjointSolver::DiscreteOrdinatesAdjointSolver(const chi::InputParameters& params)
  : lbs::DiscreteOrdinatesSolver(params)
{
  basic_options_.AddOption<std::string>("REFERENCE_RF", std::string());
}

DiscreteOrdinatesAdjointSolver::DiscreteOrdinatesAdjointSolver(const std::string& solver_name)
  : lbs::DiscreteOrdinatesSolver(solver_name)
{
  basic_options_.AddOption<std::string>("REFERENCE_RF", std::string());
}

#ifdef OPENSN_WITH_LUA
const std::vector<lbs::DiscreteOrdinatesAdjointSolver::RespFuncAndSubs>&
DiscreteOrdinatesAdjointSolver::GetResponseFunctions() const
{
  return response_functions_;
}
#endif

void
DiscreteOrdinatesAdjointSolver::Initialize()
{
  LBSSolver::Initialize();

  MakeAdjointXSs();
  InitQOIs();

  //================================================== Initialize source func
  auto src_function = std::make_shared<AdjointSourceFunction>(*this);

  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4);

  //================================================== Initialize groupsets for
  //                                                   sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  InitializeSolverSchemes(); // j
  source_event_tag_ = Chi::log.GetRepeatingEventTag("Set Source");
}

void
DiscreteOrdinatesAdjointSolver::MakeAdjointXSs()
{
  // Create adjoint cross sections
  using AdjXS = chi_physics::AdjointMGXS;

  // define the actual cross-sections
  std::map<int, XSPtr> matid_to_adj_xs_map;
  for (const auto& matid_xs_pair : matid_to_xs_map_)
  {
    const auto matid = matid_xs_pair.first;
    const auto fwd_xs = std::dynamic_pointer_cast<chi_physics::MultiGroupXS>(matid_xs_pair.second);
    matid_to_adj_xs_map[matid] = std::make_shared<AdjXS>(*fwd_xs);
  } // for each mat
  matid_to_xs_map_ = std::move(matid_to_adj_xs_map);

  // reassign transport view to adjoint cross-sections
  if (grid_ptr_->local_cells.size() == cell_transport_views_.size())
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& xs_ptr = matid_to_xs_map_[cell.material_id_];
      auto& transport_view = cell_transport_views_[cell.local_id_];

      transport_view.ReassingXS(*xs_ptr);
    }
}

void
DiscreteOrdinatesAdjointSolver::InitQOIs()
{
#ifdef OPENSN_WITH_LUA
  //============================================= Initialize QOIs
  for (auto& qoi_pair : response_functions_)
  {
    const auto& qoi_designation = qoi_pair.first;
    auto& qoi_cell_subscription = qoi_pair.second;

    for (const auto& cell : grid_ptr_->local_cells)
      if (qoi_designation.logical_volume->Inside(cell.centroid_))
        qoi_cell_subscription.push_back(cell.local_id_);

    size_t num_local_subs = qoi_cell_subscription.size();
    size_t num_globl_subs = 0;

    MPI_Allreduce(&num_local_subs, // sendbuf
                  &num_globl_subs, // recvbuf
                  1,
                  MPI_UNSIGNED_LONG_LONG, // count + datatype
                  MPI_SUM,                // operation
                  Chi::mpi.comm);         // communicator

    Chi::log.Log() << "LBAdjointSolver: Number of cells subscribed to " << qoi_designation.name
                   << " = " << num_globl_subs;
  }
#endif
}

void
DiscreteOrdinatesAdjointSolver::Execute()
{
  const std::string fname = __FUNCTION__;

  primary_ags_solver_->Setup();
  primary_ags_solver_->Solve();

  //============================================= Apply post processing
  Chi::log.Log() << "LBAdjointSolver: post-processing.";
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets_)
    for (const auto& group : groupset.groups_)
      set_group_numbers.insert(group.id_);

  const auto& m_to_ell_em_map = groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap();

  //============================================= Reorient phi-moments for
  // reverse
  //                                              angle
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_view = cell_transport_views_[cell.local_id_];
    const int num_nodes = cell_view.NumNodes();

    for (int i = 0; i < num_nodes; ++i)
    {
      for (int m = 0; m < num_moments_; ++m)
      {
        const auto& ell = m_to_ell_em_map[m].ell;

        size_t dof_map_g0 = cell_view.MapDOF(i, m, 0); // unknown map

        for (int g : set_group_numbers)
          phi_old_local_[dof_map_g0 + g] *= pow(-1.0, ell);
      } // for moment
    }   // node i
  }     // for cell

  UpdateFieldFunctions();
}

size_t
DiscreteOrdinatesAdjointSolver::AddResponseFunction(
  const std::string& qoi_name,
  std::shared_ptr<chi_mesh::LogicalVolume> logical_volume,
  const std::string& lua_function_name)
{
#ifdef OPENSN_WITH_LUA
  // Make the designation
  ResponseFunctionDesignation qoi_designation(
    qoi_name, std::move(logical_volume), lua_function_name);
  // Make empty subscriber list (will be populated during initialize)
  std::vector<size_t> cell_rf_subscriptions;

  response_functions_.emplace_back(qoi_designation, cell_rf_subscriptions);

  return response_functions_.size() - 1;
#else
  // hard exit since this is a mssing capability
  exit(-1);
#endif
}

void
DiscreteOrdinatesAdjointSolver::ExportImportanceMap(const std::string& file_name)
{
  const std::string fname = __FUNCTION__;

  //============================================= Determine cell averaged
  //                                              importance map
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets_)
    for (const auto& group : groupset.groups_)
      set_group_numbers.insert(group.id_);

  const auto& m_to_ell_em_map = groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap();

  typedef chi_math::VectorN<4> Arr4; // phi, J_x, J_y, J_z
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
        //==================================== Get multigroup p1_moments
        MGVec4 p1_moments(set_group_numbers.size(), VecDbl{0.0, 0.0, 0.0, 0.0});
        for (int m = 0; m < std::max(static_cast<int>(num_moments_), 4); ++m)
        {
          const auto& ell = m_to_ell_em_map[m].ell;
          const auto& em = m_to_ell_em_map[m].m;

          size_t dof_map_g0 = cell_view.MapDOF(i, m, 0); // unknown map

          for (int g : set_group_numbers)
          {
            if (ell == 0 and em == 0) p1_moments[g](0) = std::fabs(phi_old_local_[dof_map_g0 + g]);
            if (ell == 1 and em == 1) p1_moments[g](1) = phi_old_local_[dof_map_g0 + g];
            if (ell == 1 and em == -1) p1_moments[g](2) = phi_old_local_[dof_map_g0 + g];
            if (ell == 1 and em == 0) p1_moments[g](3) = phi_old_local_[dof_map_g0 + g];
          } // for g
        }   // for m

        nodal_p1_moments[i] = std::move(p1_moments);
      } // for node i

      //=========================================== Determine nodal average
      //                                            p1_moments
      for (int g : set_group_numbers)
      {
        chi_math::VectorN<4> cell_p1_avg(VecDbl{0.0, 0.0, 0.0, 0.0});

        double volume_total = 0.0;
        for (int i = 0; i < num_nodes; ++i)
        {
          double IntV_shapeI = fe_values.Vi_vectors[i];
          cell_p1_avg += nodal_p1_moments[i][g] * IntV_shapeI;
          volume_total += IntV_shapeI;
        } // for node i
        cell_p1_avg /= volume_total;

        cell_avg_p1_moments[cell.local_id_][g] = cell_p1_avg;
      } // for g
    }   // for cell
  }

  //============================================= Determine cell-based
  //                                              exponential-representations
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

  Chi::log.Log() << "Exporting importance map to binary file " << file_name;

  const auto locJ_io_flags = std::ofstream::binary | std::ofstream::out;
  const auto loc0_io_flags = locJ_io_flags | std::ofstream::trunc;
  const bool is_home = (Chi::mpi.location_id == 0);

  //======================================== Build header
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

  //================================================== Process each location
  uint64_t num_global_cells = grid_ptr_->GetGlobalNumberOfCells();
  for (int locationJ = 0; locationJ < Chi::mpi.process_count; ++locationJ)
  {
    Chi::log.LogAll() << "  Barrier at " << locationJ;
    Chi::mpi.Barrier();
    if (Chi::mpi.location_id != locationJ) continue;

    Chi::log.LogAll() << "  Location " << locationJ << " appending data.";

    std::ofstream file(file_name, is_home ? loc0_io_flags : locJ_io_flags);

    if (not file.is_open())
    {
      std::stringstream outstr;

      outstr << fname << ": Location " << Chi::mpi.location_id << ", failed to open file "
             << file_name;
      throw std::logic_error(outstr.str());
    }

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

  Chi::log.LogAll() << "Done exporting importance map to binary file " << file_name;
  Chi::mpi.Barrier();
}

double
DiscreteOrdinatesAdjointSolver::ComputeInnerProduct()
{
  double local_integral = 0.0;

  //============================================= Material sources
  for (const auto& cell : grid_ptr_->local_cells)
  {
    if (matid_to_src_map_.count(cell.material_id_) == 0) continue; // Skip if no src

    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& source = matid_to_src_map_[cell.material_id_];
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];

    for (const auto& group : groups_)
    {
      const int g = group.id_;
      const double Q = source->source_value_g_[g];

      if (Q > 0.0)
      {
        const int num_nodes = transport_view.NumNodes();
        for (int i = 0; i < num_nodes; ++i)
        {
          const size_t dof_map = transport_view.MapDOF(i, 0, g); // unknown map

          const double phi = phi_old_local_[dof_map];

          local_integral += Q * phi * fe_values.Vi_vectors[i];
        } // for node
      }   // check source value >0
    }     // for group
  }       // for cell

  //============================================= Point sources
  for (const auto& point_source : point_sources_)
  {
    const auto& info_list = point_source.ContainingCellsInfo();
    for (const auto& info : info_list)
    {
      const auto& cell = grid_ptr_->local_cells[info.cell_local_id];
      const auto& transport_view = cell_transport_views_[cell.local_id_];
      const auto& source_strength = point_source.Strength();
      const auto& shape_values = info.shape_values;

      for (const auto& group : groups_)
      {
        const int g = group.id_;
        const double S = source_strength[g] * info.volume_weight;

        if (S > 0.0)
        {
          const int num_nodes = transport_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            const size_t dof_map = transport_view.MapDOF(i, 0, g); // unknown
                                                                   // map

            const double phi_i = phi_old_local_[dof_map];

            local_integral += S * phi_i * shape_values[i];
          } // for node
        }   // check source value >0
      }     // for group
    }       // for cell
  }         // for point source

  double global_integral = 0.0;

  MPI_Allreduce(&local_integral,  // sendbuf
                &global_integral, // recvbuf
                1,
                MPI_DOUBLE,     // count, datatype
                MPI_SUM,        // op
                Chi::mpi.comm); // comm

  return global_integral;
}

} // namespace lbs
