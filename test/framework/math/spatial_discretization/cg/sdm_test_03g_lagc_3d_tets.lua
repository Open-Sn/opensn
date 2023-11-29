dofile("mesh_3d_tets.lua")

unit_tests.math_SDM_Test01_Continuous
({
  sdm_type = "LagrangeC",
  --export_vtk = true
});
