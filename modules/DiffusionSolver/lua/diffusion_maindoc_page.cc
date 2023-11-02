/** \defgroup LuaDiffusion Diffusion
 * \ingroup LuaModules
 *
 *
 * Please consult the whitepaper for this solver (<a
 * href="whitepages/DiffusionSolver/DiffusionSolver.pdf">
 * Diffusion Whitepaper</a>). This solver solves the general diffusion equation
 * of the form
 *
 * \f{equation}{
 * -\nabla D \nabla \phi + \sigma_a \phi = q.
 * \f}
 *
 * Given a discretization method this solver will assemble the matrix \f$ A \f$
 * using
 * a system of linear equations. The solver can operate on 1D slabs, 2D polygon
 * meshes and 3D polyhedron meshes. Currently two spatial discretization schemes
 * are supported, Piecewise Linear Continous (PWLC) and Piecewise Linear
 * Discontinous (PWLD) using the Modified Interior Penalty (MIP) method. The
 * solver is fully parallel (using PETSc).
 *
 * Boundary conditions are specified by referencing unique boundary id's. i.e.
 * assuming one uses the extruder mesher the template boundary would've been
 * assigned to index 0 and this will be transfused to all the boundaries in the
 * lateral periphery. The top and bottom boundaries are flat for extruded
 * geometries and therefore are assigned their own id's. In the extruded case
 * the last boundary index is always assigned to the top boundary and
 * second-to-last boundary index is always assigned to the bottom boundary.
 * Boundary types are specified using the ChiLua::chiDiffusionSetProperty
 * function call, using the BOUNDARY_TYPE property index.
 *
 * The materials and source values are for now obtained from materials
 * associated with cells. The default property mapping is shown below. If
 * the property at index [1] (which is mapped to the source value q) is not
 * available the source will default to a constant value of 1.0.
 *
 *  \image html "DiffusionMatProp.png" width=500px
 *
 * To change the mapping of the properties the user needs to make a call to
 * ChiLua::chiDiffusionSetProperty using the PROPERTY_D_MAP, PROPERTY_Q_MAP or
 * PROPERTY_SIGMAA_MAP property index.
 *
 * By default the solver populates a scalar field function which is the solution
 * \f$ \phi \f$.
 *
 * */