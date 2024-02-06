--############################################### Create cross sections
xs = {}

for m=0,6 do
    xs[tostring(m)] = PhysicsTransportXSCreate()
end

-- GMesh mesh
PhysicsTransportXSSet(xs["0"],OPENSN_XSFILE,"materials/XS_water.xs")
PhysicsTransportXSSet(xs["1"],OPENSN_XSFILE,"materials/XS_UO2.xs")
PhysicsTransportXSSet(xs["2"],OPENSN_XSFILE,"materials/XS_7pMOX.xs")
PhysicsTransportXSSet(xs["3"],OPENSN_XSFILE,"materials/XS_guide_tube.xs")
PhysicsTransportXSSet(xs["4"],OPENSN_XSFILE,"materials/XS_4_3pMOX.xs")
PhysicsTransportXSSet(xs["5"],OPENSN_XSFILE,"materials/XS_8_7pMOX.xs")
PhysicsTransportXSSet(xs["6"],OPENSN_XSFILE,"materials/XS_fission_chamber.xs")

water_xs = PhysicsTransportXSGet(xs["0"])

num_groups = water_xs["num_groups"]
Log(LOG_0,"Num groups: "..tostring(num_groups))

--############################################### Create materials
materials = {}
for m=0,6 do
    key = tostring(m)
    materials[key] = PhysicsAddMaterial("Material_"..key)
    PhysicsMaterialAddProperty(key,TRANSPORT_XSECTIONS)
    PhysicsMaterialSetProperty(key,TRANSPORT_XSECTIONS, EXISTING,xs[key])
end
