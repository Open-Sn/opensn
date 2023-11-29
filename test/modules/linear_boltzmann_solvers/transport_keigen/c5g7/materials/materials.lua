--############################################### Create cross sections
xs = {}

for m=0,6 do
    xs[tostring(m)] = chiPhysicsTransportXSCreate()
end

-- GMesh mesh
chiPhysicsTransportXSSet(xs["0"],CHI_XSFILE,"materials/XS_water.cxs")
chiPhysicsTransportXSSet(xs["1"],CHI_XSFILE,"materials/XS_UO2.cxs")
chiPhysicsTransportXSSet(xs["2"],CHI_XSFILE,"materials/XS_7pMOX.cxs")
chiPhysicsTransportXSSet(xs["3"],CHI_XSFILE,"materials/XS_guide_tube.cxs")
chiPhysicsTransportXSSet(xs["4"],CHI_XSFILE,"materials/XS_4_3pMOX.cxs")
chiPhysicsTransportXSSet(xs["5"],CHI_XSFILE,"materials/XS_8_7pMOX.cxs")
chiPhysicsTransportXSSet(xs["6"],CHI_XSFILE,"materials/XS_fission_chamber.cxs")

water_xs = chiPhysicsTransportXSGet(xs["0"])

num_groups = water_xs["num_groups"]
chiLog(LOG_0,"Num groups: "..tostring(num_groups))

--############################################### Create materials
materials = {}
for m=0,6 do
    key = tostring(m)
    materials[key] = PhysicsAddMaterial("Material_"..key)
    PhysicsMaterialAddProperty(key,TRANSPORT_XSECTIONS)
    PhysicsMaterialSetProperty(key,TRANSPORT_XSECTIONS, EXISTING,xs[key])
end
