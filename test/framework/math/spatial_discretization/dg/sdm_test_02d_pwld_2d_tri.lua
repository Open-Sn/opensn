dofile("mesh_2d_tri.lua")

unit_tests.math_SDM_Test02_DisContinuous({
  mesh = grid,
  sdm_type = "PWLD",
})
