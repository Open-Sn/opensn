lv1 = logvol.RPPLogicalVolume.Create({ xmin = -12.0, xmax = 12.0 })

print("lv1 test 1:", lv1:Inside({ x = -13.0, y = 0.1, z = 0.1 }))
print("lv1 test 2:", lv1:Inside({ x = -11.0, y = 0.1, z = 0.1 }))
print("lv1 test 3:", lv1:Inside({ x = -11.0, y = -1.0, z = -1.0 }))

lv2 = logvol.RPPLogicalVolume.Create({ xmin = -12.0, xmax = 12.0, infy = True, infz = True })
print()
print("lv2 test 1:", lv2:Inside({ x = -13.0, y = 0.1, z = 0.1 }))
print("lv2 test 2:", lv2:Inside({ x = -11.0, y = 0.1, z = 0.1 }))
print("lv2 test 3:", lv2:Inside({ x = -11.0, y = -1.0, z = -1.0 }))
