lv1 = chi_mesh.RPPLogicalVolume.Create({xmin = -12.0, xmax = 12.0})

print("lv1 test 1:", chiLogicalVolumePointSense(lv1, -13.0, 0.1, 0.1))
print("lv1 test 2:", chiLogicalVolumePointSense(lv1, -11.0, 0.1, 0.1))
print("lv1 test 3:", chiLogicalVolumePointSense(lv1, -11.0, -1.0, -1.0))

lv2 = chi_mesh.RPPLogicalVolume.Create({xmin = -12.0, xmax = 12.0,
                                        infy = true, infz = true})
print()
print("lv2 test 1:", chiLogicalVolumePointSense(lv2, -13.0, 0.1, 0.1))
print("lv2 test 2:", chiLogicalVolumePointSense(lv2, -11.0, 0.1, 0.1))
print("lv2 test 3:", chiLogicalVolumePointSense(lv2, -11.0, -1.0, -1.0))