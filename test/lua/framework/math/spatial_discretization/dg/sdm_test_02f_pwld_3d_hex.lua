dofile("mesh_3d_hex.lua")

unit_tests.math_SDM_Test02_DisContinuous({
  mesh = grid,
  sdm_type = "PWLD",
  --export_vtk = true,
})
