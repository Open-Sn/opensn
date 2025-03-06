dofile("mesh_2d_ortho.lua")

unit_tests.math_SDM_Test02_DisContinuous({
  mesh = grid,
  sdm_type = "PWLD",
})
