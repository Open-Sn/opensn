dofile("mesh_3d_prism.lua")

unit_tests.math_SDM_Test01_Continuous({
  mesh = grid,
  sdm_type = "PWLC",
})
