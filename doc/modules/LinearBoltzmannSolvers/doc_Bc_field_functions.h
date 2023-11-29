/** \page LBSFieldFunctions Field Functions

LBS solvers create a number of field functions. Here is the logic for it:
- Field functions are by default named `phi_gXXX_mYYY` where
  `XXX` is a zero padded 3 digit number for the group
  number. `YYY` is a zero padded 3 digit number for the moment number. Numbers
  spanning beyond the 3 digits will have no zero padding and will be
  represented normally (the 3 digit padding just helps in visualization cases).
Example:
Suppose this is a 3D simulation, 2 groups, scattering order of 1
(resulting in 4 moments)
\code
phys1 = LBSCreateSolver()
SolverAddRegion(phys1,region1)
--
-- Add Groupset construction here
--
LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
LBSSetProperty(phys1,SCATTERING_ORDER,1)
--
chiLBSInitialize(phys1)
chiLBSExecute(phys1)
--
fflist,count = LBSGetScalarFieldFunctionList(phys1)
\endcode

will create field functions
\code
phi_g000_m000
phi_g000_m001
phi_g000_m002
phi_g000_m003
phi_g001_m000
phi_g001_m001
phi_g001_m002
phi_g001_m003
\endcode

We can get the scalar field functions with a call to
`LBSGetScalarFieldFunctionList` which will return a table with the field
function handles of only the scalar fields. E.g.,
\code
fflist = LBSGetScalarFieldFunctionList(phys1)
\endcode
with `fflist` containing handles to
\code
phi_g000_m000
phi_g001_m000
\endcode

Additionally, LBS can create a power generation field function with the
default name `power_generation`.
\code
chiLBSSetOptions(phys1,
{
  spatial_discretization = "pwld",
  scattering_order = 2,
  power_field_function_on = true,
  power_default_kappa = 1.0,
  power_normalization = -1.0, --Negative means its disabled
})
\endcode

All of the field function can be supplied with a prefix, either using the
solver name
\code
chiLBSSetOptions(phys1,
{
  spatial_discretization = "pwld",
  scattering_order = 2,
  power_field_function_on = true,
  power_default_kappa = 1.0,
  power_normalization = -1.0, --Negative means its disabled

  field_function_prefix_option = "solver_name"
})
\endcode

or by setting a different prefix
\code
chiLBSSetOptions(phys1,
{
  spatial_discretization = "pwld",
  scattering_order = 2,
  power_field_function_on = true,
  power_default_kappa = 1.0,
  power_normalization = -1.0, --Negative means its disabled

  field_function_prefix = "prefixx"
})
\endcode

\ingroup LBSUtilities*/
