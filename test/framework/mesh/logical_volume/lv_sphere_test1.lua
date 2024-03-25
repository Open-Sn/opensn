lv1 = logvol.SphereLogicalVolume.Create({r = 1.3, x=1.0, y=-1.0, z=2.0})

print("lv1 test 1:", logvol.PointSense(lv1, 1.0,-1.0,2.0))
print("lv1 test 2:", logvol.PointSense(lv1, 1.0+1.3,-1.0,2.0))
print("lv1 test 3:", logvol.PointSense(lv1, 1.0+1.301,-1.0,2.0))
