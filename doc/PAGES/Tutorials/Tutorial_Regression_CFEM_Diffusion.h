/** \page Tutorials_CFEM_Diffusion  CFEM Diffusion Test-cases as Tutorials

The following CFEM diffusion test-cases are built on the same 2D square geometry:
1. **[CFEM_Diffusion_2D_1a_linear](https://github.com/chi-tech/chi-tech/blob/development/tests/CFEM_Diffusion/CFEM_Diffusion_2D_1a_linear.lua)**
   - *description:* pure scatterer, no volumetric source, boundary source, solution is linear
   - *feature of potential interest:*
      - material properties a Lua functions of position and material ID,
      - use of chiLogicalVolume to define boundaries,
      - field function and export of FF line to python.
2. **[CFEM_Diffusion_2D_2a_DirichletBCs](https://github.com/chi-tech/chi-tech/blob/development/tests/CFEM_Diffusion/CFEM_Diffusion_2D_2a_DirichletBCs.lua)**
   - *description:* diffusion-reaction on a 2-region domain, with Dirichlet BCs
   - *feature of potential interest:*
         - use of chiLogicalVolume to define material regions,
         - material properties a Lua functions of position and material ID.
3. **[CFEM_Diffusion_2D_3b_analytical_coef2](https://github.com/chi-tech/chi-tech/blob/development/tests/CFEM_Diffusion/CFEM_Diffusion_2D_3b_analytical_coef2.lua)**
   - *description:* diffusion-reaction problem with space-dependent material properties
   - *feature of potential interest:*
         - use of chiLogicalVolume to define material regions,
         - material properties a Lua functions of position and material ID.

 */