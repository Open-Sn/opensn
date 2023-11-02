#include "lbs_solver.h"

#include "ChiObjectFactory.h"

namespace lbs
{

// ##################################################################
RegisterSyntaxBlock(/*namespace_name=*/lbs,
                    /*block_name=*/OptionsBlock,
                    /*syntax_function=*/LBSSolver::OptionsBlock);

chi::InputParameters
LBSSolver::OptionsBlock()
{
  chi::InputParameters params;

  params.SetGeneralDescription("Set options from a large list of parameters");
  params.SetDocGroup("LBSUtilities");

  // clang-format off
  params.AddOptionalParameter("spatial_discretization", "pwld",
    "What spatial discretization to use. Currently only `\"pwld\"` "
    "is supported");
  params.AddOptionalParameter("scattering_order",1,
  "Defines the level of harmonic expansion for the scattering source.");
  params.AddOptionalParameter("sweep_eager_limit",32'000,
  "The eager limit to be used in message size during sweep initialization.\n"
  " This expects to be followed by a size in bytes (Max 64,0000)See note below."
  "\\n\\n"
  " ###Note on the Eager limit\n"
  "The eager limit is the message size limit before which non-blocking MPI send"
  "calls will execute without waiting for a matching receive call. The limit is"
  "platform dependent but in general 64 kb. Some systems have 32 kb as a limit"
  "and therefore we use that as a default limit in ChiTech. There is a fine"
  "interplay between message size and the shear amount of messages that will be"
  "sent. In general smaller messages tend to be more efficient, however, when"
  "there are too many small messages being sent around the communication system"
  "on the given platform will start to suffer. One can gain a small amount of"
  "parallel efficiency by lowering this limit, however, there is a point where"
  "the parallel efficiency will actually get worse so use with caution.");
  params.AddOptionalParameter("read_restart_data",false,
  "Flag indicating whether restart data is to be read.");
  params.AddOptionalParameter("read_restart_folder_name","YRestart",
  "Folder name to use when reading restart data.");
  params.AddOptionalParameter("read_restart_file_base","restart",
  "File base name to use when reading restart data.");
  params.AddOptionalParameter("write_restart_data",false,
  "Flag indicating whether restart data is to be written.");
  params.AddOptionalParameter("write_restart_folder_name","YRestart",
  "Folder name to use when writing restart data.");
  params.AddOptionalParameter("write_restart_file_base","restart",
  "File base name to use when writing restart data.");
  params.AddOptionalParameter("write_restart_interval",30.0,
  "Interval at which restart data is to be written. Currently not implemented.");
  params.AddOptionalParameter("use_precursors",false,
  "Flag for using delayed neutron precursors.");
  params.AddOptionalParameter("use_source_moments",false,
  "Flag for ignoring fixed sources and selectively using source moments "
  "obtained elsewhere.");
  params.AddOptionalParameter("save_angular_flux",false,
  "Flag indicating whether angular fluxes are to be stored or not.");
  params.AddOptionalParameter("verbose_inner_iterations",true,
  "Flag to control verbosity of inner iterations.");
  params.AddOptionalParameter("verbose_outer_iterations",true,
  "Flag to control verbosity of across-groupset iterations.");
  params.AddOptionalParameter("verbose_ags_iterations",false,
  "Flag to control verbosity of across-groupset iterations.");
  params.AddOptionalParameter("power_field_function_on",false,
  "Flag to control the creation of the power generation field function. If set "
  "to `true` then a field function will be created with the general name "
  "`<solver_name>_power_generation`.");
  params.AddOptionalParameter("power_default_kappa",3.20435e-11,
  "Default `kappa` value (Energy released per fission) to use for power "
  "generation when cross sections do not have `kappa` values. Default: "
  "3.20435e-11 Joule (corresponding to 200 MeV per fission).");
  params.AddOptionalParameter("power_normalization",-1.0,
  "Power normalization factor to use. Supply a negative or zero number to turn "
  "this off.");
  params.AddOptionalParameter("field_function_prefix_option","prefix",
  "Prefix option on field function names. Default: `\"prefix\"`. Can be "
  "`\"prefix\"` or "
  "`\"solver_name\"`. "
  "By default this option is `\"prefix\"` which means it uses the designated "
  "\"prefix\" (another option), however, that is defaulted to nothing. "
  "Therefore, default behavior is to export flux moment fields functions as "
  "`phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit group number and "
  "similarly for `YYY`.");
  params.AddOptionalParameter("field_function_prefix","",
  "Prefix to use on all field functions. Default: `\"\"`. "
  "By default this option is empty but if specified then flux moments will "
  "exported as `prefix_phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit "
  "group number and similarly for `YYY`. The underscore after \"prefix\" is "
  "added automatically.");
  params.AddOptionalParameterArray("boundary_conditions",
  {},
  "A table contain sub-tables for each boundary specification.");
  params.LinkParameterToBlock("boundary_conditions",
                              "lbs::BoundaryOptionsBlock");

  using namespace chi_data_types;
  params.ConstrainParameterRange("spatial_discretization",
      AllowableRangeList::New({"pwld"}));

  params.ConstrainParameterRange("field_function_prefix_option",
    AllowableRangeList::New({"prefix", "solver_name"}));
  // clang-format on

  return params;
}

// ##################################################################
RegisterSyntaxBlock(/*namespace_in_lua=*/lbs,
                    /*name_in_lua=*/BoundaryOptionsBlock,
                    /*syntax_function=*/LBSSolver::BoundaryOptionsBlock);

chi::InputParameters
LBSSolver::BoundaryOptionsBlock()
{
  chi::InputParameters params;

  // clang-format off
  params.SetGeneralDescription(
    "Set options for boundary conditions. See \\ref LBSBCs");
  params.SetDocGroup("LBSUtilities");

  params.AddRequiredParameter<std::string>("name",
  "Boundary name that identifies the specific boundary");
  params.AddRequiredParameter<std::string>("type",
  "Boundary type specification.");

  params.AddOptionalParameterArray<double>("group_strength", {},
  "Required only if `type` is `\"incident_isotropic\"`. An array of isotropic "
  "strength per group");

  params.AddOptionalParameter("function_name", "",
  "Text name of the lua function to be called for this boundary condition. For"
  " more on this boundary condition type.");

  using namespace chi_data_types;
  params.ConstrainParameterRange("name", AllowableRangeList::New({
  "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}));

  params.ConstrainParameterRange("type", AllowableRangeList::New({
  "vacuum", "incident_isotropic", "reflecting",
  "incident_anisotropic_heterogeneous"}));
  // clang-format on

  return params;
}

// ##################################################################
void
LBSSolver::SetOptions(const chi::InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();

  for (size_t p = 0; p < user_params.NumParameters(); ++p)
  {
    const auto& spec = user_params.GetParam(p);

    if (spec.Name() == "spatial_discretization")
    {
      auto sdm_name = spec.GetValue<std::string>();
      if (sdm_name == "pwld")
        Options().sd_type = chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
    }

    else if (spec.Name() == "scattering_order")
      Options().scattering_order = spec.GetValue<int>();

    else if (spec.Name() == "sweep_eager_limit")
      Options().sweep_eager_limit = spec.GetValue<int>();

    else if (spec.Name() == "read_restart_data")
      Options().read_restart_data = spec.GetValue<bool>();

    else if (spec.Name() == "read_restart_folder_name")
      Options().read_restart_folder_name = spec.GetValue<std::string>();

    else if (spec.Name() == "read_restart_file_base")
      Options().read_restart_file_base = spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_data")
      Options().write_restart_data = spec.GetValue<bool>();

    else if (spec.Name() == "write_restart_folder_name")
      Options().write_restart_folder_name = spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_file_base")
      Options().write_restart_file_base = spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_interval")
      Options().write_restart_interval = spec.GetValue<double>();

    else if (spec.Name() == "use_precursors")
      Options().use_precursors = spec.GetValue<bool>();

    else if (spec.Name() == "use_source_moments")
      Options().use_src_moments = spec.GetValue<bool>();

    else if (spec.Name() == "save_angular_flux")
      Options().save_angular_flux = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_inner_iterations")
      Options().verbose_inner_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_ags_iterations")
      Options().verbose_ags_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_outer_iterations")
      Options().verbose_outer_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "power_field_function_on")
      Options().power_field_function_on = spec.GetValue<bool>();

    else if (spec.Name() == "power_default_kappa")
      Options().power_default_kappa = spec.GetValue<double>();

    else if (spec.Name() == "power_normalization")
      Options().power_normalization = spec.GetValue<double>();

    else if (spec.Name() == "field_function_prefix_option")
    {
      Options().field_function_prefix_option = spec.GetValue<std::string>();
    }

    else if (spec.Name() == "field_function_prefix")
      Options().field_function_prefix = spec.GetValue<std::string>();

    else if (spec.Name() == "boundary_conditions")
    {
      spec.RequireBlockTypeIs(chi::ParameterBlockType::ARRAY);

      for (size_t b = 0; b < spec.NumParameters(); ++b)
      {
        auto bndry_params = BoundaryOptionsBlock();
        bndry_params.AssignParameters(spec.GetParam(b));

        SetBoundaryOptions(bndry_params);
      }
    }
  } // for p
}

// ##################################################################
void
LBSSolver::SetBoundaryOptions(const chi::InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const auto& user_params = params.ParametersAtAssignment();
  const auto boundary_name = user_params.GetParamValue<std::string>("name");
  const auto bndry_type = user_params.GetParamValue<std::string>("type");

  const std::map<std::string, uint64_t> supported_bndry_names = {
    {"xmin", 1}, {"xmax", 0}, {"ymin", 3}, {"ymax", 2}, {"zmin", 5}, {"zmax", 4}};
  const auto bid = supported_bndry_names.at(boundary_name);
  const std::map<std::string, lbs::BoundaryType> type_list = {
    {"vacuum", BoundaryType::VACUUM},
    {"incident_isotropic", BoundaryType::INCIDENT_ISOTROPIC},
    {"reflecting", BoundaryType::REFLECTING},
    {"incident_anisotropic_heterogeneous", BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS}};

  const auto type = type_list.at(bndry_type);
  switch (type)
  {
    case BoundaryType::VACUUM:
    case BoundaryType::REFLECTING:
    {
      BoundaryPreferences()[bid] = {type};
      break;
    }
    case BoundaryType::INCIDENT_ISOTROPIC:
    {
      if (not user_params.Has("group_strength"))
      {
        std::string message = fname;
        message += ":boundary_conditions:"
                   "type=\"incident_isotropic\" requires parameter "
                   "\"group_strength\".";

        throw std::invalid_argument(message);
      }
      user_params.RequireParameterBlockTypeIs("group_strength", chi::ParameterBlockType::ARRAY);

      const auto group_strength = user_params.GetParamVectorValue<double>("group_strength");
      BoundaryPreferences()[bid] = {type, group_strength};
      break;
    }
    case BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS:
    {
      if (not user_params.Has("function_name"))
      {
        std::string message = fname;
        message += ":boundary_conditions:"
                   "type=\"incident_anisotropic_heterogeneous\" requires "
                   "parameter \"function_name\".";

        throw std::invalid_argument(message);
      }
      const auto bndry_function_name = user_params.GetParamValue<std::string>("function_name");

      BoundaryPreferences()[bid] = {type, {}, bndry_function_name};
      break;
    }
  }
}

} // namespace lbs