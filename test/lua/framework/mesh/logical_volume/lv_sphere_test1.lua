lv1 = logvol.SphereLogicalVolume.Create({ r = 1.3, x = 1.0, y = -1.0, z = 2.0 })

print("lv1 test 1:", lv1:Inside({ x = 1.0, y = -1.0, z = 2.0 }))
print("lv1 test 2:", lv1:Inside({ x = 1.0 + 1.3, y = -1.0, z = 2.0 }))
print("lv1 test 3:", lv1:Inside({ x = 1.0 + 1.301, y = -1.0, z = 2.0 }))
