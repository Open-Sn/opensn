--[[ @doc
# The Rest of the Transport Simulation Input

## Materials and Sources
We create two materials and add two properties to it:
+ TRANSPORT_XSECTIONS for the transport cross sections, and
+ ISOTROPIC_MG_SOURCE for the isotropic volumetric source

We assign cross sections per material by loading the file containing the cross sections.

We create lua tables containing the volumetric multigroup sources and assign them to their
 respective material by passing an array.
--]]
-- Add materials
materials = {}
materials[1] = PhysicsAddMaterial("material-1");
materials[2] = PhysicsAddMaterial("material-2");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  OPENSN_XSFILE,"xs_1g_MatA.xs")
PhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
  OPENSN_XSFILE,"xs_1g_MatB.xs")

PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
PhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
num_groups = 1
src1={}
src2={}
for g=1,num_groups do
  src1[g] = 0.0
  src2[g] = 300.0
end
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src1)
PhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src2)

--[[ @doc
## Angular Quadrature

We call a product Gauss-Legendre-Chebyshev quadrature and pass the number of **positive** polar cosines
(here ```npolar = 2```) and the number of azimuthal subdivisions in **one quadrant** (```nazimu = 1```).
This creates a 3D angular quadrature.
--]]
-- Setup the Angular Quadrature
nazimu = 4
npolar = 2
pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,nazimu,npolar)

-- Setup LBS parameters
--[[ @doc
## lbs_solver_rest Linear Boltzmann Solver
### Options for the Linear Boltzmann Solver (LBS)
In the LBS block, we provide
+ the number of energy groups,
+ the groupsets (with 0-indexing), the handle for the angular quadrature, the angle aggregation, the solver type, tolerances, and other solver options.
--]]
lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 0},
      angular_quadrature_handle = pquad,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    }
  }
}
--[[ @doc
### Further Options for the Linear Boltzmann Solver
In the LBS options, we pass the maximum scattering order to be employed (should be less than the one supplied the cross section file)
--]]
lbs_options =
{
  scattering_order = 0,
}
--[[ @doc
### Putting the Linear Boltzmann Solver Together
We create the physics solver, initialize it, and execute it,
--]]
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys})

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

