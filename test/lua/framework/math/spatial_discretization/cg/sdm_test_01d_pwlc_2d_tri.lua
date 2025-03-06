dofile("mesh_2d_tri.lua")

unit_tests.math_SDM_Test01_Continuous({
  mesh = grid,
  sdm_type = "PWLC",
})
