#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/lbs_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/source_functions/source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/diffusion_dfem_solver/iterative_methods/mip_wgs_context2.h"
#include "framework/object_factory.h"

namespace opensn
{
namespace lbs
{

OpenSnRegisterObjectInNamespace(lbs, DiffusionDFEMSolver);

InputParameters
DiffusionDFEMSolver::GetInputParameters()
{
  InputParameters params = LBSSolver::GetInputParameters();

  params.SetClassName("DiffusionDFEMSolver");
  params.SetDocGroup("lbs__LBSSolver");

  params.ChangeExistingParamToOptional("name", "LBSDiffusionDFEMSolver");

  return params;
}

DiffusionDFEMSolver::DiffusionDFEMSolver(const InputParameters& params) : LBSSolver(params)
{
}

DiffusionDFEMSolver::~DiffusionDFEMSolver()
{
  for (auto& groupset : groupsets_)
    CleanUpTGDSA(groupset);
}

void
DiffusionDFEMSolver::Initialize()
{
  options_.scattering_order = 0; // overwrite any setting otherwise
  LBSSolver::Initialize();

  auto src_function = std::make_shared<SourceFunction>(*this);

  // Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4, _5);

  // Initialize groupsets preconditioning
  for (auto& groupset : groupsets_)
    InitTGDSA(groupset);

  LBSSolver::InitializeSolverSchemes();
}

void
DiffusionDFEMSolver::InitializeWGSSolvers()
{
  // Initialize groupset solvers
  gs_mip_solvers_.assign(groupsets_.size(), nullptr);
  const size_t num_groupsets = groupsets_.size();
  for (size_t gs = 0; gs < num_groupsets; ++gs)
  {
    const auto& groupset = groupsets_[gs];

    // Make UnknownManager
    const size_t gs_G = groupset.groups_.size();
    opensn::UnknownManager uk_man;
    uk_man.AddUnknown(UnknownType::VECTOR_N, gs_G);

    // Make boundary conditions
    std::map<uint64_t, BoundaryCondition> bcs;
    for (auto& [bid, lbs_bndry] : sweep_boundaries_)
    {
      if (lbs_bndry->Type() == BoundaryType::REFLECTING)
        bcs[bid] = {BCType::ROBIN, {0.0, 1.0, 0.0}};
      else if (lbs_bndry->Type() == BoundaryType::ISOTROPIC)
      {
        const bool has_bndry_preference = boundary_preferences_.count(bid) > 0;
        if (not has_bndry_preference)
          bcs[bid] = {BCType::ROBIN, {0.25, 0.5}};

        const auto& bpref = boundary_preferences_.at(bid);
        const bool is_vaccuum = bpref.type == BoundaryType::VACUUM;
        if (is_vaccuum)
          bcs[bid] = {BCType::ROBIN, {0.25, 0.5}};
        else
          throw std::logic_error("Dirichlet boundary conditions not supported"
                                 "for diffusion solvers.");
      }
    } // for sweep-boundary

    // Make xs map
    typedef Multigroup_D_and_sigR MGXS;
    typedef std::map<int, Multigroup_D_and_sigR> MatID2XSMap;
    MatID2XSMap matid_2_mgxs_map;
    for (const auto& matid_xs_pair : matid_to_xs_map_)
    {
      const auto& mat_id = matid_xs_pair.first;
      const auto& xs = matid_xs_pair.second;

      const auto& diffusion_coeff = xs->DiffusionCoefficient();
      const auto& sigma_r = xs->SigmaRemoval();

      std::vector<double> Dg(gs_G, 0.0);
      std::vector<double> sigR(gs_G, 0.0);

      size_t g = 0;
      for (size_t gprime = groupset.groups_.front().id_; gprime <= groupset.groups_.back().id_;
           ++gprime)
      {
        Dg[g] = diffusion_coeff[gprime];
        sigR[g] = sigma_r[gprime];
        ++g;
      } // for g

      matid_2_mgxs_map.insert(std::make_pair(mat_id, MGXS{Dg, sigR}));
    }

    // Create solver
    const auto& sdm = *discretization_;

    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(TextName() + "_WGSolver"),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       unit_cell_matrices_,
                                                       true); // verbosity

    solver->options.residual_tolerance = groupset.wgdsa_tol_;
    solver->options.max_iters = groupset.wgdsa_max_iters_;
    solver->options.verbose = groupset.wgdsa_verbose_;
    solver->options.additional_options_string = groupset.wgdsa_string_;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    gs_mip_solvers_[gs] = solver;
  } // for groupset

  wgs_solvers_.clear(); // this is required
  for (auto& groupset : groupsets_)
  {

    auto mip_wgs_context_ptr = std::make_shared<MIPWGSContext2>(
      *this,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES | SUPPRESS_WG_SCATTER,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options_.verbose_inner_iterations);

    auto wgs_solver = std::make_shared<WGSLinearSolver>(mip_wgs_context_ptr);

    wgs_solvers_.push_back(wgs_solver);
  } // for groupset
}

} // namespace lbs
} // namespace opensn
