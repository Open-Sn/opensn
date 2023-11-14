dofile("mesh_3d_prism.lua")

chi_unit_tests.chi_math_SDM_Test02_DisContinuous
({
  sdm_type = "PWLD",
  --export_vtk = true,
});
