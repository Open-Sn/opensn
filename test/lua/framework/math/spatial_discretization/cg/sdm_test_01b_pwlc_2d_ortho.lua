dofile("mesh_2d_ortho.lua")

unit_tests.math_SDM_Test01_Continuous({
  mesh = grid,
  sdm_type = "PWLC",
})
