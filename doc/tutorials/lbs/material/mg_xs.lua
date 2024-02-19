--[[ @doc
# Multigroup Cross Sections

## Background

OpenSn is not provided with cross-section libraries. Users are expected to supply their own multigroup cross-section data.
One may use other open-source software to generate this data (e.g., NJOY, Dragon, OpenMC).

The cross sections are read from a plain text file. The ```OPENSN_XSFILE``` format  of that file is as follows:
```
# Add descriptive comments as appropriate
NUM_GROUPS ng
NUM_MOMENTS nmom

SIGMA_T_BEGIN
0 value
.
.
ng-1 value
SIGMA_T_END

SIGMA_A_BEGIN
0 value
.
.
ng-1 value
SIGMA_A_END

TRANSFER_MOMENTS_BEGIN
# Add descriptive comments as appropriate
M_GPRIME_G_VAL 0 0 0 value
.
M_GPRIME_G_VAL moment gprime g value
.
M_GPRIME_G_VAL nmom-1 ng-1 ng-1 value
TRANSFER_MOMENTS_END

```

## Mesh

A simple orthogonal 2D mesh.

We create a right parallelepiped logical volume that contains the entire mesh and we assign a 0 for material ID to 
all cells found inside the logical volume.

--]]
-- Setup the mesh
nodes={}
n_cells=4
length=2.
xmin = -length/2.
dx = length/n_cells
for i=1,(n_cells+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create
({
    node_sets = {nodes,nodes},
})

mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)

--[[ @doc

## Materials

We create a material and add two properties to it:
+ TRANSPORT_XSECTIONS for the transport cross sections, and
+ ISOTROPIC_MG_SOURCE for the isotropic volumetric source
--]]
-- Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Material_A");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

--[[ @doc

## Cross Sections

We assign the cross sections to the material by loading the file containing the cross sections.
--]]
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        OPENSN_XSFILE,"xs_1g_MatA.xs")

--[[ @doc

## Volumetric Source

We create a lua table containing the volumetric multigroup source and assign it to the material by passing that array.
--]]
num_groups = 1
src={}
for g=1,num_groups do
    src[g] = 1.0
end
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--[[ @doc

## Angular Quadrature

We call a product Gauss-Legendre-Chebyshev quadrature and pass the number of **positive** polar cosines
(here ```npolar = 2```) and the number of azimuthal subdivisions in **one quadrant** (```nazimu = 1```).
This creates a 3D angular quadrature.

We finish by optimizing the quadrature to only use the positive hemisphere for 2D simulations.
--]]
-- Setup the Angular Quadrature
nazimu = 1
npolar = 2
pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,nazimu,npolar)
OptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

--[[ @doc

## Linear Boltzmann Solver

In the LBS block, we provide
+ the number of energy groups,
+ the groupsets (with 0-indexing), the handle for the angular quadrature, the angle aggregation, the solver type,
tolerances, and other solver options.

In the LBS options, we pass the maximum scattering order to be employed (should be less than the one supplied the
cross-section file)

We then create the physics solver, initialize it, and execute it.

--]]
-- Setup LBS parameters
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

lbs_options =
{
    scattering_order = 0,
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys})

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

--[[ @doc

## Post-Processing via Field Functions

We extract the scalar flux (i.e., the first entry in the field function list; recall that lua indexing starts at 1)
and export it to a VTK file whose name is supplied by the user.

--]]
-- Get field functions
fflist,count = LBSGetScalarFieldFunctionList(phys)
vtk_basename = "mg_xs"
ExportFieldFunctionToVTK(fflist[1],vtk_basename)
