dofile("mesh_3d_hex.lua")

chi_unit_tests.math_SDM_Test01_Continuous
({
  sdm_type = "LagrangeC",
  --export_vtk = true
});
