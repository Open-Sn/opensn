lv1 = logvol.RPPLogicalVolume.Create({xmin = -12.0, xmax = 12.0})

print("lv1 test 1:", logvol.PointSense(lv1, {x = -13.0, y = 0.1, z = 0.1}))
print("lv1 test 2:", logvol.PointSense(lv1, {x = -11.0, y = 0.1, z = 0.1}))
print("lv1 test 3:", logvol.PointSense(lv1, {x = -11.0, y = -1.0, z = -1.0}))

lv2 = logvol.RPPLogicalVolume.Create({xmin = -12.0, xmax = 12.0,
                                        infy = true, infz = true})
print()
print("lv2 test 1:", logvol.PointSense(lv2, {x = -13.0, y = 0.1, z = 0.1}))
print("lv2 test 2:", logvol.PointSense(lv2, {x = -11.0, y = 0.1, z = 0.1}))
print("lv2 test 3:", logvol.PointSense(lv2, {x = -11.0, y = -1.0, z = -1.0}))
