#include "opensn/framework/chi_lua.h"

#include "opensn/modules/MGDiffusion/mg_diffusion_solver.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"

namespace mg_diffusion::mgd_lua_utils
{

//#############################################################################
/** Sets a property of a Diffusion solver. Please also consult the whitepaper
 * for the Diffusion solver (<a
 * href="../../whitepages/DiffusionSolver/DiffusionSolver.pdf">
 * Diffusion Whitepaper</a>)

\n\n Additional basic options can be set as indicated in \ref
LuaDiffusionBasicOptions

\param SolverHandle int Handle to an existing diffusion solver.
\param PropertyName string Name for a specific property.
\param Values varying Number of inputs associated with the index.<br>

##_

###PropertyName\n
"boundary_type"\n
 Boundary type. Expects boundary index then <B>BoundaryTypeName</B>
 then type value.\n\n

\code
chiCFEMMGDiffusionSetBCProperty(solver,"boundary_type",bdID,"vacuum")
\endcode

### BoundaryTypeName
reflecting\n
 Reflecting boundary conditions. Synonymous with Neumann with a
 derivative of 0.0.
             \f[ -D \hat{n}\cdot \nabla \phi = 0 \f]\n\n

neumann\n
 Constant derivative boundary condition. Expects to be followed
 by a constant \f$ f \f$ representing
                    \f[ -D \hat{n}\cdot \nabla \phi = f \f]\n\n

vacuum\n
 Vacuum boundary conditions. More appropriate to neutron diffusion.
   \f[ \frac{1}{4}\phi + \frac{1}{2} D \hat{n}\cdot \nabla \phi = 0 \f]\n\n

robin\n
 Robin boundary condition of the form
                   \f[ a \phi + b D \hat{n}\cdot \nabla \phi = f \f]\n\n

\ingroup LuaDiffusion
\author Jean*/
int
chiCFEMMGDiffusionSetBCProperty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, num_args, 2);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  //==========================sss=================== Get solver
  LuaCheckNumberValue(fname, L, 1);
  const int solver_index = lua_tonumber(L, 1);

  auto& solver = Chi::GetStackItem<mg_diffusion::Solver>(Chi::object_stack, solver_index, fname);

  //============================================= Get property index
  LuaCheckStringValue(fname, L, 2);
  const std::string property_name = lua_tostring(L, 2);

  //============================================= Handle properties
  if (property_name == "boundary_type")
  {
    if (num_args < 4)
    {
      Chi::log.Log0Error() << "Invalid amount of arguments used in"
                           << " chiCFEMMGDiffusionsetBCproperty(...,\"boundary_type\".... "
                           << " At least 4 arguments are expected.";
      Chi::Exit(EXIT_FAILURE);
    }
    LuaCheckNumberValue(fname, L, 3);
    const int bound_index = lua_tonumber(L, 3);

    LuaCheckStringValue(fname, L, 4);
    const std::string type_name = lua_tostring(L, 4);

    if (type_name == "reflecting") // ------------- REFLECTING
    {
      if (num_args != 4)
      {
        Chi::log.Log0Error() << "Invalid amount of arguments used in"
                             << " chiCFEMMGDiffusionsetBCproperty(...,\"boundary_type\","
                             << bound_index << ",\"reflecting\". "
                             << " 4 arguments are expected.";
        Chi::Exit(EXIT_FAILURE);
      }

      mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = mg_diffusion::BoundaryType::Reflecting;

      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      Chi::log.Log() << "Boundary " << bound_index << " set as "
                     << "Reflecting.";
    }
    else if (type_name == "dirichlet") // ------------- DIRICHLET
    {
      Chi::log.Log0Error() << "Dirichlet BC is not supported in multigroup diffusion "
                           << "(chiCFEMMGDiffusionSetBCProperty).";
      Chi::Exit(EXIT_FAILURE);
    }
    else if (type_name == "neumann") // ------------- NEUMANN
    {
      if (num_args != 5)
      {
        Chi::log.Log0Error() << "Invalid amount of arguments used in"
                             << " chiCFEMMGDiffusionsetBCproperty(...,\"boundary_type\","
                             << bound_index << ",\"neumann\". "
                             << " 5 arguments are expected.";
        Chi::Exit(EXIT_FAILURE);
      }
      // check lua tables
      LuaCheckTableValue(fname, L, 5);
      std::vector<double> f_values;
      LuaPopulateVectorFrom1DArray(fname, L, 5, f_values);
      // add the other multigroup vectors to finish the BC
      unsigned int ng = f_values.size();
      std::vector<double> a_values(ng, 0.0);
      std::vector<double> b_values(ng, 1.0);

      mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = mg_diffusion::BoundaryType::Neumann;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      Chi::log.Log() << "Boundary " << bound_index << " set as "
                     << "Neumann with D_g grad(u_g) dot n = f_g";
    }
    else if (type_name == "vacuum") // ------------- VACUUM
    {
      if (num_args != 4)
      {
        Chi::log.Log0Error() << "Invalid amount of arguments used in"
                             << " chiCFEMMGDiffusionsetBCproperty(...,\"boundary_type\","
                             << bound_index << ",\"vacuum\". "
                             << " 4 arguments are expected.";
        Chi::Exit(EXIT_FAILURE);
      }
      // dummy-sized values until we now num_group later, after solver init
      std::vector<double> a_values(1, 0.25);
      std::vector<double> b_values(1, 0.5);
      std::vector<double> f_values(1, 0.0);

      mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = mg_diffusion::BoundaryType::Vacuum;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      Chi::log.Log() << "Boundary " << bound_index << " set as "
                     << "Vacuum.";
    }
    else if (type_name == "robin") // ------------- ROBIN
    {
      if (num_args != 7)
      {
        Chi::log.Log0Error() << "Invalid amount of arguments used in"
                             << " chiCFEMMGDiffusionSetBCProperty(...,\"boundary_type\","
                             << bound_index << ",\"robin\". "
                             << " 7 arguments are expected.";
        Chi::Exit(EXIT_FAILURE);
      }
      // check lua tables
      LuaCheckTableValue(fname, L, 5);
      LuaCheckTableValue(fname, L, 6);
      LuaCheckTableValue(fname, L, 7);
      std::vector<double> a_values, b_values, f_values;
      LuaPopulateVectorFrom1DArray(fname, L, 5, a_values);
      LuaPopulateVectorFrom1DArray(fname, L, 6, b_values);
      LuaPopulateVectorFrom1DArray(fname, L, 7, f_values);

      mg_diffusion::Solver::BoundaryInfo bndry_info;
      bndry_info.first = mg_diffusion::BoundaryType::Robin;
      bndry_info.second = {a_values, b_values, f_values};
      solver.boundary_preferences_.insert(std::make_pair(bound_index, bndry_info));

      Chi::log.Log() << "Boundary " << bound_index << " set as Robin";
    }
    else
    {
      Chi::log.LogAllError() << "Unsupported boundary type encountered in call to "
                             << "chiCFEMMGDiffusionSetBCProperty(..,\"boundary_type\",.. :"
                             << type_name;
      Chi::Exit(EXIT_FAILURE);
    }
  }
  else // wrong property_name
  {
    Chi::log.Log0Error() << "Invalid property in chiCFEMMGDiffusionSetBCProperty.";
    Chi::Exit(EXIT_FAILURE);
  }
  return 0;
} // end of chiCFEMMGDiffusionSetBCProperty

} // namespace mg_diffusion::mgd_lua_utils
