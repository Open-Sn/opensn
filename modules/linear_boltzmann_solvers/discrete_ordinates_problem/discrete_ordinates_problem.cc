// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/vacuum_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/isotropic_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aah_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/tgdsa.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/logging/log.h"
#include "framework/logging/log_exceptions.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <iomanip>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, DiscreteOrdinatesProblem);

DiscreteOrdinatesProblem::DiscreteOrdinatesProblem(const std::string& name,
                                                   std::shared_ptr<MeshContinuum> grid_ptr)
  : LBSProblem(name, grid_ptr)
{
}

InputParameters
DiscreteOrdinatesProblem::GetInputParameters()
{
  InputParameters params = LBSProblem::GetInputParameters();

  params.SetClassName("DiscreteOrdinatesProblem");

  params.ChangeExistingParamToOptional("name", "LBSDiscreteOrdinatesProblem");

  params.AddOptionalParameterArray(
    "directions_sweep_order_to_print",
    std::vector<int>(),
    "List of direction id's for which sweep ordering info is to be printed.");

  params.AddOptionalParameter(
    "sweep_type", "AAH", "The sweep type to use for sweep operatorations.");

  params.ConstrainParameterRange("sweep_type", AllowableRangeList::New({"AAH", "CBC"}));

  return params;
}

std::shared_ptr<DiscreteOrdinatesProblem>
DiscreteOrdinatesProblem::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<DiscreteOrdinatesProblem>("lbs::DiscreteOrdinatesProblem", params);
}

DiscreteOrdinatesProblem::DiscreteOrdinatesProblem(const InputParameters& params)
  : LBSProblem(params),
    verbose_sweep_angles_(params.GetParamVectorValue<int>("directions_sweep_order_to_print")),
    sweep_type_(params.GetParamValue<std::string>("sweep_type"))
{
  if (use_gpus_ && sweep_type_ == "CBC")
  {
    log.Log0Warning() << "Sweep computation on GPUs has not yet been supported for CBC. "
                      << "Falling back to CPU sweep.\n";
    use_gpus_ = false;
  }
}

DiscreteOrdinatesProblem::~DiscreteOrdinatesProblem()
{
  CALI_CXX_MARK_FUNCTION;

  for (auto& groupset : groupsets_)
  {
    WGDSA::CleanUp(groupset);
    TGDSA::CleanUp(groupset);

    // Reset sweep orderings
    if (groupset.angle_agg != nullptr)
      groupset.angle_agg->angle_set_groups.clear();
  }
}

std::pair<size_t, size_t>
DiscreteOrdinatesProblem::GetNumPhiIterativeUnknowns()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::GetNumPhiIterativeUnknowns");
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_global_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  size_t num_local_psi_dofs = 0;
  size_t num_global_psi_dofs = 0;
  for (auto& groupset : groupsets_)
  {
    const auto num_delayed_psi_info = groupset.angle_agg->GetNumDelayedAngularDOFs();
    num_local_psi_dofs += num_delayed_psi_info.first;
    num_global_psi_dofs += num_delayed_psi_info.second;
  }

  const size_t num_local_dofs = num_local_phi_dofs + num_local_psi_dofs;
  const size_t num_global_dofs = num_global_phi_dofs + num_global_psi_dofs;

  return {num_local_dofs, num_global_dofs};
}

const std::map<uint64_t, std::shared_ptr<SweepBoundary>>&
DiscreteOrdinatesProblem::GetSweepBoundaries() const
{
  return sweep_boundaries_;
}

void
DiscreteOrdinatesProblem::Initialize()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::Initialize");

  LBSProblem::Initialize();

  // Make face histogram
  grid_face_histogram_ = grid_->MakeGridFaceHistogram();

  const auto grid_dim = grid_->GetDimension();
  for (auto& groupset : groupsets_)
  {
    const auto quad_dim = groupset.quadrature->GetDimension();
    OpenSnInvalidArgumentIf(grid_dim != quad_dim,
                            "Dimensionality of quadrature set (" + std::to_string(quad_dim) +
                              ") for groupset #" + std::to_string(groupset.id) +
                              " does not match dimensionality of mesh (" +
                              std::to_string(grid_dim) + ").");
  }

  // Initialize source func
  using namespace std::placeholders;
  auto src_function = std::make_shared<SourceFunction>(*this);
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4); // NOLINT

  // Initialize groupsets for sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    WGDSA::Init(*this, groupset);
    TGDSA::Init(*this, groupset);
  }
  InitializeSolverSchemes();
}

void
DiscreteOrdinatesProblem::InitializeBoundaries()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeBoundaries");

  // Determine boundary-ids involved in the problem
  std::set<uint64_t> global_unique_bids_set;
  {
    std::set<uint64_t> local_unique_bids_set;
    for (const auto& cell : grid_->local_cells)
      for (const auto& face : cell.faces)
        if (not face.has_neighbor)
          local_unique_bids_set.insert(face.neighbor_id);

    std::vector<uint64_t> local_unique_bids(local_unique_bids_set.begin(),
                                            local_unique_bids_set.end());
    std::vector<uint64_t> recvbuf;
    mpi_comm.all_gather(local_unique_bids, recvbuf);

    global_unique_bids_set = local_unique_bids_set; // give it a head start

    for (uint64_t bid : recvbuf)
      global_unique_bids_set.insert(bid);
  }

  // Initialize default incident boundary
  const size_t G = num_groups_;

  sweep_boundaries_.clear();
  for (uint64_t bid : global_unique_bids_set)
  {
    const bool has_no_preference = boundary_preferences_.count(bid) == 0;
    const bool has_not_been_set = sweep_boundaries_.count(bid) == 0;
    if (has_no_preference and has_not_been_set)
    {
      sweep_boundaries_[bid] = std::make_shared<VacuumBoundary>(G);
    } // defaulted
    else if (has_not_been_set)
    {
      const auto& bndry_pref = boundary_preferences_.at(bid);
      const auto& mg_q = bndry_pref.isotropic_mg_source;

      if (bndry_pref.type == LBSBoundaryType::VACUUM)
        sweep_boundaries_[bid] = std::make_shared<VacuumBoundary>(G);
      else if (bndry_pref.type == LBSBoundaryType::ISOTROPIC)
        sweep_boundaries_[bid] = std::make_shared<IsotropicBoundary>(G, mg_q);
      else if (bndry_pref.type == LBSBoundaryType::REFLECTING)
      {
        // Locally check all faces, that subscribe to this boundary,
        // have the same normal
        const double EPSILON = 1.0e-12;
        std::unique_ptr<Vector3> n_ptr = nullptr;
        for (const auto& cell : grid_->local_cells)
          for (const auto& face : cell.faces)
            if (not face.has_neighbor and face.neighbor_id == bid)
            {
              if (not n_ptr)
                n_ptr = std::make_unique<Vector3>(face.normal);
              if (std::fabs(face.normal.Dot(*n_ptr) - 1.0) > EPSILON)
                throw std::logic_error(
                  "LBSProblem: Not all face normals are, within tolerance, locally the "
                  "same for the reflecting boundary condition requested");
            }

        // Now check globally
        const int local_has_bid = n_ptr != nullptr ? 1 : 0;
        const Vector3 local_normal = local_has_bid ? *n_ptr : Vector3(0.0, 0.0, 0.0);

        std::vector<int> locJ_has_bid(opensn::mpi_comm.size(), 1);
        std::vector<double> locJ_n_val(opensn::mpi_comm.size() * 3, 0.0);

        mpi_comm.all_gather(local_has_bid, locJ_has_bid);
        std::vector<double> lnv = {local_normal.x, local_normal.y, local_normal.z};
        mpi_comm.all_gather(lnv.data(), 3, locJ_n_val.data(), 3);

        Vector3 global_normal;
        for (int j = 0; j < opensn::mpi_comm.size(); ++j)
        {
          if (locJ_has_bid[j])
          {
            int offset = 3 * j;
            const double* n = &locJ_n_val[offset];
            const Vector3 locJ_normal(n[0], n[1], n[2]);

            if (local_has_bid)
              if (std::fabs(local_normal.Dot(locJ_normal) - 1.0) > EPSILON)
                throw std::logic_error(
                  "LBSProblem: Not all face normals are, within tolerance, globally the "
                  "same for the reflecting boundary condition requested");

            global_normal = locJ_normal;
          }
        }

        sweep_boundaries_[bid] = std::make_shared<ReflectingBoundary>(
          G, global_normal, MapGeometryTypeToCoordSys(options_.geometry_type));
      }
    } // non-defaulted
  } // for bndry id
}

void
DiscreteOrdinatesProblem::InitializeWGSSolvers()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeWGSSolvers");

  // Determine max size and number of matrices along sweep front
  for (auto& groupset : groupsets_)
  {
    // Max groupset size
    max_groupset_size_ = std::max(max_groupset_size_, groupset.groups.size());

    for (auto& angle_set_group : groupset.angle_agg->angle_set_groups)
    {
      for (auto& angleset : angle_set_group.GetAngleSets())
      {
        // Max level size
        const auto& spds = angleset->GetSPDS();
        const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
        for (const auto& level : levelized_spls)
          max_level_size_ = std::max(max_level_size_, level.size());

        // Max angleset size
        max_angleset_size_ = std::max(max_angleset_size_, angleset->GetAngleIndices().size());
      }
    }
  }

  wgs_solvers_.clear(); // this is required
  for (auto& groupset : groupsets_)
  {
    std::shared_ptr<SweepChunk> sweep_chunk = SetSweepChunk(groupset);
    auto sweep_wgs_context_ptr = std::make_shared<SweepWGSContext>(
      *this,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options_.verbose_inner_iterations,
      sweep_chunk);

    if (groupset.iterative_method == LinearSystemSolver::IterativeMethod::CLASSIC_RICHARDSON)
    {
      wgs_solvers_.push_back(std::make_shared<ClassicRichardson>(
        sweep_wgs_context_ptr, options_.verbose_inner_iterations));
    }
    else
      wgs_solvers_.push_back(std::make_shared<WGSLinearSolver>(sweep_wgs_context_ptr));
  }
}

void
DiscreteOrdinatesProblem::ReorientAdjointSolution()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::ReorientAdjointSolution");

  for (const auto& groupset : groupsets_)
  {
    int gs = groupset.id;

    // Moment map for flux moments
    const auto& moment_map = groupset.quadrature->GetMomentToHarmonicsIndexMap();

    // Angular flux info
    auto& psi = psi_new_local_[gs];
    const auto& uk_man = groupset.psi_uk_man_;

    // Build reversed angle mapping
    std::map<int, int> reversed_angle_map;
    if (options_.save_angular_flux)
    {
      const auto& omegas = groupset.quadrature->omegas;
      const auto num_gs_angles = omegas.size();

      // Go through angles until all are paired
      std::set<size_t> visited;
      for (size_t idir = 0; idir < num_gs_angles; ++idir)
      {
        // Skip if already encountered
        if (visited.count(idir) > 0)
          continue;

        bool found = true;
        for (size_t jdir = 0; jdir < num_gs_angles; ++jdir)
        {
          // Angles are opposite if their sum is zero
          const auto sum = grid_->GetDimension() == 1
                             ? Vector3(0.0, 0.0, omegas[idir].z + omegas[jdir].z)
                             : omegas[idir] + omegas[jdir];
          const bool opposite = sum.NormSquare() < 1.0e-8;

          // Add opposites to mapping
          if (opposite)
          {
            found = true;
            reversed_angle_map[idir] = jdir;

            visited.insert(idir);
            visited.insert(jdir);
          }
        } // for angle n

        OpenSnLogicalErrorIf(not found,
                             "Opposing angle for " + omegas[idir].PrintStr() + " in groupset " +
                               std::to_string(gs) + " not found.");

      } // for angle m
    } // if saving angular flux

    const auto num_gs_groups = groupset.groups.size();
    const auto gsg_i = groupset.groups.front().id;
    const auto gsg_f = groupset.groups.back().id;

    for (const auto& cell : grid_->local_cells)
    {
      const auto& transport_view = cell_transport_views_[cell.local_id];
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        // Reorient flux moments
        //
        // Because flux moments are integrated angular fluxes, the
        // angular flux and spherical harmonics must be evaluated at
        // opposite angles in the quadrature integration. Taking advantage
        // of the even/odd nature of the spherical harmonics, i.e.
        // Y_{\ell,m}(-\Omega) = (-1)^\ell Y_{\ell,m}(\Omega), the flux
        // moments must be multiplied by (-1)^\ell.
        for (size_t imom = 0; imom < num_moments_; ++imom)
        {
          const auto& ell = moment_map[imom].ell;
          const auto dof_map = transport_view.MapDOF(i, imom, 0);

          for (int g = gsg_i; g <= gsg_f; ++g)
          {
            phi_new_local_[dof_map + g] *= std::pow(-1.0, ell);
            phi_old_local_[dof_map + g] *= std::pow(-1.0, ell);
          } // for group g
        } // for moment m

        // Reorient angular flux
        if (options_.save_angular_flux)
        {
          for (const auto& [idir, jdir] : reversed_angle_map)
          {
            const auto dof_map =
              std::make_pair(discretization_->MapDOFLocal(cell, i, uk_man, idir, 0),
                             discretization_->MapDOFLocal(cell, i, uk_man, jdir, 0));

            for (size_t gsg = 0; gsg < num_gs_groups; ++gsg)
              std::swap(psi[dof_map.first + gsg], psi[dof_map.second + gsg]);
          }
        }
      } // for node i
    } // for cell

  } // for groupset
}

void
DiscreteOrdinatesProblem::ZeroOutflowBalanceVars(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::ZeroOutflowBalanceVars");

  for (const auto& cell : grid_->local_cells)
    for (int f = 0; f < cell.faces.size(); ++f)
      for (auto& group : groupset.groups)
        cell_transport_views_[cell.local_id].ZeroOutflow(f, group.id);
}

void
DiscreteOrdinatesProblem::InitializeSweepDataStructures()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeSweepDataStructures");

  log.Log() << program_timer.GetTimeString() << " Initializing sweep datastructures.\n";

  // Define sweep ordering groups
  quadrature_unq_so_grouping_map_.clear();
  std::map<std::shared_ptr<AngularQuadrature>, bool> quadrature_allow_cycles_map_;
  for (auto& groupset : groupsets_)
  {
    if (quadrature_unq_so_grouping_map_.count(groupset.quadrature) == 0)
    {
      quadrature_unq_so_grouping_map_[groupset.quadrature] = AssociateSOsAndDirections(
        grid_, *groupset.quadrature, groupset.angleagg_method, options_.geometry_type);
    }

    if (quadrature_allow_cycles_map_.count(groupset.quadrature) == 0)
      quadrature_allow_cycles_map_[groupset.quadrature] = groupset.allow_cycles;
  }

  // Build sweep orderings
  quadrature_spds_map_.clear();
  if (sweep_type_ == "AAH")
  {
    // Creating an AAH SPDS can be an expensive operation. We break it up into multiple phases so
    // so that we can distribute the work across MPI ranks:
    // 1) Initialize the SPDS for each angleset. This is done by all ranks.
    // 2) For each SPDS, generate the feedback arc set (FAS) for the global sweep graph. Angelsets
    //    are distributed as evenly as possible across MPI ranks.
    // 3) Gather the FAS for each SPDS on all ranks.
    // 4) Build the global sweep task dependency graph (TDG) for each SPDS.

    // Initalize SPDS. All ranks initialize a SPDS for each angleset.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Initializing AAH SPDS.";
    for (const auto& [quadrature, info] : quadrature_unq_so_grouping_map_)
    {
      int id = 0;
      const auto& unique_so_groupings = info.first;
      for (const auto& so_grouping : unique_so_groupings)
      {
        if (so_grouping.empty())
          continue;

        const size_t master_dir_id = so_grouping.front();
        const auto& omega = quadrature->omegas[master_dir_id];
        const auto new_swp_order = std::make_shared<AAH_SPDS>(
          id, omega, this->grid_, quadrature_allow_cycles_map_[quadrature]);
        quadrature_spds_map_[quadrature].push_back(new_swp_order);
        ++id;
      }
    }

    // Generate the global sweep FAS for each SPDS. This is an expensive operation. It is
    // distributed via MPI so that multiple MPI ranks can compute the FAS for one or more SPDS
    // independently.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep FAS for each SPDS.";
    for (const auto& quadrature : quadrature_spds_map_)
    {
      for (const auto& spds : quadrature.second)
      {
        auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
        auto id = aah_spds->GetId();
        if (opensn::mpi_comm.rank() == (id % opensn::mpi_comm.size()))
          aah_spds->BuildGlobalSweepFAS();
      }
    }

    // Communicate the FAS for each SPDS to all ranks.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Gather FAS for each SPDS.";
    std::vector<int> local_edges_to_remove;
    for (const auto& quadrature : quadrature_spds_map_)
    {
      for (const auto& spds : quadrature.second)
      {
        auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
        auto id = aah_spds->GetId();
        if ((id % opensn::mpi_comm.size()) == opensn::mpi_comm.rank())
        {
          auto edges_to_remove = aah_spds->GetGlobalSweepFAS();
          local_edges_to_remove.push_back(id);
          local_edges_to_remove.push_back(static_cast<int>(edges_to_remove.size()));
          local_edges_to_remove.insert(
            local_edges_to_remove.end(), edges_to_remove.begin(), edges_to_remove.end());
        }
      }
    }

    int local_size = static_cast<int>(local_edges_to_remove.size());
    std::vector<int> receive_counts(opensn::mpi_comm.size(), 0);
    std::vector<int> displacements(opensn::mpi_comm.size(), 0);
    mpi_comm.all_gather(local_size, receive_counts);

    int total_size = 0;
    for (auto i = 0; i < receive_counts.size(); ++i)
    {
      displacements[i] = total_size;
      total_size += receive_counts[i];
    }

    std::vector<int> global_edges_to_remove(total_size, 0);
    mpi_comm.all_gather(
      local_edges_to_remove, global_edges_to_remove, receive_counts, displacements);

    // Unpack the gathered data and update SPDS on all ranks.
    int offset = 0;
    while (offset < global_edges_to_remove.size())
    {
      int spds_id = global_edges_to_remove[offset++];
      int num_edges = global_edges_to_remove[offset++];
      std::vector<int> edges;
      edges.reserve(num_edges);
      for (auto i = 0; i < num_edges; ++i)
        edges.emplace_back(global_edges_to_remove[offset++]);

      for (const auto& quadrature : quadrature_spds_map_)
      {
        for (const auto& spds : quadrature.second)
        {
          auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
          if (aah_spds->GetId() == spds_id)
          {
            aah_spds->SetGlobalSweepFAS(edges);
            break;
          }
        }
      }
    }

    // Build TDG for each SPDS on all ranks.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep TDGs.";
    for (const auto& quadrature : quadrature_spds_map_)
      for (const auto& spds : quadrature.second)
        std::static_pointer_cast<AAH_SPDS>(spds)->BuildGlobalSweepTDG();

    // Print ghosted sweep graph if requested
    if (not verbose_sweep_angles_.empty())
    {
      for (const auto& quadrature : quadrature_spds_map_)
      {
        for (const auto& spds : quadrature.second)
        {
          for (const int dir_id : verbose_sweep_angles_)
          {
            auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
            if (aah_spds->GetId() == dir_id)
              aah_spds->PrintGhostedGraph();
          }
        }
      }
    }
  }
  else if (sweep_type_ == "CBC")
  {
    // Build SPDS
    for (const auto& [quadrature, info] : quadrature_unq_so_grouping_map_)
    {
      const auto& unique_so_groupings = info.first;
      for (const auto& so_grouping : unique_so_groupings)
      {
        if (so_grouping.empty())
          continue;

        const size_t master_dir_id = so_grouping.front();
        const auto& omega = quadrature->omegas[master_dir_id];
        const auto new_swp_order =
          std::make_shared<CBC_SPDS>(omega, this->grid_, quadrature_allow_cycles_map_[quadrature]);
        quadrature_spds_map_[quadrature].push_back(new_swp_order);
      }
    }
  }
  else
    OpenSnInvalidArgument("Unsupported sweep type \"" + sweep_type_ + "\"");

  opensn::mpi_comm.barrier();

  // Build FLUDS templates
  quadrature_fluds_commondata_map_.clear();
  if (sweep_type_ == "AAH")
  {
    for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
    {
      for (const auto& spds : spds_list)
      {
        quadrature_fluds_commondata_map_[quadrature].push_back(
          std::make_unique<AAH_FLUDSCommonData>(
            grid_nodal_mappings_, *spds, *grid_face_histogram_));
      }
    }
  }
  else if (sweep_type_ == "CBC")
  {
    for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
    {
      for (const auto& spds : spds_list)
      {
        quadrature_fluds_commondata_map_[quadrature].push_back(
          std::make_unique<CBC_FLUDSCommonData>(*spds, grid_nodal_mappings_));
      }
    }
  }

  log.Log() << program_timer.GetTimeString() << " Done initializing sweep datastructures.\n";
}

std::pair<UniqueSOGroupings, DirIDToSOMap>
DiscreteOrdinatesProblem::AssociateSOsAndDirections(const std::shared_ptr<MeshContinuum> grid,
                                                    const AngularQuadrature& quadrature,
                                                    const AngleAggregationType agg_type,
                                                    const GeometryType lbs_geo_type)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::AssociateSOsAndDirections");

  const std::string fname = __FUNCTION__;

  // Checks
  if (quadrature.omegas.empty())
    throw std::logic_error(fname + ": Quadrature with no omegas cannot be used.");
  if (quadrature.weights.empty())
    throw std::logic_error(fname + ": Quadrature with no weights cannot be used.");

  // Build groupings
  UniqueSOGroupings unq_so_grps;
  switch (agg_type)
  {
    // Single
    // The easiest aggregation type. Every direction
    // either has/is assumed to have a unique sweep
    // ordering. Hence there is only group holding ALL
    // the direction indices.
    case AngleAggregationType::SINGLE:
    {
      const size_t num_dirs = quadrature.omegas.size();
      for (size_t n = 0; n < num_dirs; ++n)
        unq_so_grps.push_back({n});
      break;
    } // case agg_type SINGLE

      // Polar
      // The following conditions allow for polar
      // angle aggregation.
    case AngleAggregationType::POLAR:
    {
      // Check geometry types
      if (grid->GetType() != ORTHOGONAL and grid->GetDimension() != 2 and not grid->Extruded())
        throw std::logic_error(
          fname + ": The simulation is using polar angle aggregation for which only certain "
                  "geometry types are supported, i.e., ORTHOGONAL, 2D or 3D EXTRUDED.");

      // Check quadrature type
      const auto quad_type = quadrature.GetType();
      if (quad_type != AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(fname + ": The simulation is using polar angle aggregation for "
                                       "which only Product-type quadratures are supported.");

      // Process Product Quadrature
      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        const auto num_azi = product_quad.azimu_ang.size();
        const auto num_pol = product_quad.polar_ang.size();

        // Make two separate list of polar angles
        // One upward-pointing and one downward
        std::vector<size_t> upward_polar_ids;
        std::vector<size_t> dnward_polar_ids;
        for (size_t p = 0; p < num_pol; ++p)
          if (product_quad.polar_ang[p] > M_PI_2)
            upward_polar_ids.push_back(p);
          else
            dnward_polar_ids.push_back(p);

        // Define lambda working for both upward and dnward polar-ids
        /**Lambda to convert indices and push it onto unq_so_grps.*/
        auto MapPolarAndAzimuthalIDs =
          [&product_quad, &unq_so_grps](const DirIDs& polar_ids, const size_t azimuthal_id)
        {
          DirIDs dir_ids;
          dir_ids.reserve(polar_ids.size());
          for (const size_t p : polar_ids)
            dir_ids.push_back(product_quad.GetAngleNum(p, azimuthal_id));
          unq_so_grps.push_back(std::move(dir_ids));
        };

        // Stack id's for all azimuthal angles
        for (size_t a = 0; a < num_azi; ++a)
        {
          if (not upward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(upward_polar_ids, a);
          if (not dnward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(dnward_polar_ids, a);
        } // for azi-id a

      } // try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(
          fname + ": Casting the angular quadrature to the product quadrature base, failed.");
      }

      break;
    } // case agg_type POLAR

      // Azimuthal
    case AngleAggregationType::AZIMUTHAL:
    {
      // Check geometry types
      if (lbs_geo_type != GeometryType::ONED_SPHERICAL and
          lbs_geo_type != GeometryType::TWOD_CYLINDRICAL)
        throw std::logic_error(
          fname + ": The simulation is using azimuthal angle aggregation for which only "
                  "ONED_SPHERICAL or TWOD_CYLINDRICAL derived geometry types are supported.");

      // Check quadrature type
      const auto quad_type = quadrature.GetType();
      if (quad_type != AngularQuadratureType::ProductQuadrature)
        throw std::logic_error(fname + ": The simulation is using azimuthal angle aggregation for "
                                       "which only Product-type quadratures are supported.");

      // Process Product Quadrature
      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        for (const auto& dir_set : product_quad.GetDirectionMap())
        {
          std::vector<unsigned int> group1;
          std::vector<unsigned int> group2;
          for (const auto& dir_id : dir_set.second)
            if (quadrature.abscissae[dir_id].phi > M_PI_2)
              group1.push_back(dir_id);
            else
              group2.push_back(dir_id);

          DirIDs group1_ids(group1.begin(), group1.end());
          DirIDs group2_ids(group2.begin(), group2.end());

          unq_so_grps.push_back(std::move(group1_ids));
          unq_so_grps.push_back(std::move(group2_ids));
        }
      } // try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(
          fname + ": Casting the angular quadrature to the product quadrature base, failed.");
      }

      break;
    }
    default:
      throw std::invalid_argument(fname + ": Called with UNDEFINED angle "
                                          "aggregation type.");
  } // switch angle aggregation type

  // Map directions to sweep orderings
  DirIDToSOMap dir_id_to_so_map;
  {
    size_t so_grouping_id = 0;
    for (const auto& so_grouping : unq_so_grps)
    {
      for (const size_t dir_id : so_grouping)
        dir_id_to_so_map[dir_id] = so_grouping_id;

      ++so_grouping_id;
    } // for so_grouping
  } // map scope

  return {unq_so_grps, dir_id_to_so_map};
}

void
DiscreteOrdinatesProblem::InitFluxDataStructures(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitFluxDataStructures");

  const auto& quadrature_sweep_info = quadrature_unq_so_grouping_map_[groupset.quadrature];

  const auto& unique_so_groupings = quadrature_sweep_info.first;
  const auto& dir_id_to_so_map = quadrature_sweep_info.second;

  const size_t gs_num_grps = groupset.groups.size();

  // Passing the sweep boundaries to the angle aggregation
  groupset.angle_agg =
    std::make_shared<AngleAggregation>(sweep_boundaries_, num_groups_, groupset.quadrature, grid_);

  AngleSetGroup angle_set_group;
  size_t angle_set_id = 0;
  for (const auto& so_grouping : unique_so_groupings)
  {
    const size_t master_dir_id = so_grouping.front();
    const size_t so_id = dir_id_to_so_map.at(master_dir_id);

    const auto& sweep_ordering = quadrature_spds_map_[groupset.quadrature][so_id];
    const auto& fluds_common_data = *quadrature_fluds_commondata_map_[groupset.quadrature][so_id];

    // Compute direction subsets
    const auto dir_subsets = MakeSubSets(so_grouping.size(), groupset.master_num_ang_subsets);

    for (const auto& dir_ss_info : dir_subsets)
    {
      const auto& dir_ss_begin = dir_ss_info.ss_begin;
      const auto& dir_ss_end = dir_ss_info.ss_end;
      const auto& dir_ss_size = dir_ss_info.ss_size;

      std::vector<size_t> angle_indices(dir_ss_size, 0);
      {
        size_t k = 0;
        for (size_t n = dir_ss_begin; n <= dir_ss_end; ++n)
          angle_indices[k++] = so_grouping[n];
      }

      if (sweep_type_ == "AAH")
      {
        std::shared_ptr<FLUDS> fluds =
          std::make_shared<AAH_FLUDS>(gs_num_grps,
                                      angle_indices.size(),
                                      dynamic_cast<const AAH_FLUDSCommonData&>(fluds_common_data));

        auto angle_set = std::make_shared<AAH_AngleSet>(angle_set_id++,
                                                        gs_num_grps,
                                                        *sweep_ordering,
                                                        fluds,
                                                        angle_indices,
                                                        sweep_boundaries_,
                                                        options_.max_mpi_message_size,
                                                        *grid_local_comm_set_,
                                                        use_gpus_);

        angle_set_group.GetAngleSets().push_back(angle_set);
      }
      else if (sweep_type_ == "CBC")
      {
        OpenSnLogicalErrorIf(not options_.save_angular_flux,
                             "When using sweep_type \"CBC\" then "
                             "\"save_angular_flux\" must be true.");
        std::shared_ptr<FLUDS> fluds =
          std::make_shared<CBC_FLUDS>(gs_num_grps,
                                      angle_indices.size(),
                                      dynamic_cast<const CBC_FLUDSCommonData&>(fluds_common_data),
                                      psi_new_local_[groupset.id],
                                      groupset.psi_uk_man_,
                                      *discretization_);

        auto angle_set = std::make_shared<CBC_AngleSet>(angle_set_id++,
                                                        gs_num_grps,
                                                        *sweep_ordering,
                                                        fluds,
                                                        angle_indices,
                                                        sweep_boundaries_,
                                                        *grid_local_comm_set_,
                                                        use_gpus_);

        angle_set_group.GetAngleSets().push_back(angle_set);
      }
      else
        OpenSnInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
    } // for an_ss
  } // for so_grouping

  groupset.angle_agg->angle_set_groups.push_back(std::move(angle_set_group));

  if (options_.verbose_inner_iterations)
    log.Log() << program_timer.GetTimeString() << " Initialized angle aggregation.";

  opensn::mpi_comm.barrier();
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::SetSweepChunk(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::SetSweepChunk");

  if (sweep_type_ == "AAH")
  {
    auto sweep_chunk = std::make_shared<AAHSweepChunk>(grid_,
                                                       *discretization_,
                                                       unit_cell_matrices_,
                                                       cell_transport_views_,
                                                       densities_local_,
                                                       phi_new_local_,
                                                       psi_new_local_[groupset.id],
                                                       q_moments_local_,
                                                       groupset,
                                                       block_id_to_xs_map_,
                                                       num_moments_,
                                                       max_cell_dof_count_,
                                                       *this,
                                                       max_level_size_,
                                                       max_groupset_size_,
                                                       max_angleset_size_,
                                                       use_gpus_);

    return sweep_chunk;
  }
  else if (sweep_type_ == "CBC")
  {
    auto sweep_chunk = std::make_shared<CBCSweepChunk>(phi_new_local_,
                                                       psi_new_local_[groupset.id],
                                                       grid_,
                                                       *discretization_,
                                                       unit_cell_matrices_,
                                                       cell_transport_views_,
                                                       densities_local_,
                                                       q_moments_local_,
                                                       groupset,
                                                       block_id_to_xs_map_,
                                                       num_moments_,
                                                       max_cell_dof_count_);

    return sweep_chunk;
  }
  else
    OpenSnLogicalError("Unsupported sweep_type_ \"" + sweep_type_ + "\"");
}

} // namespace opensn
