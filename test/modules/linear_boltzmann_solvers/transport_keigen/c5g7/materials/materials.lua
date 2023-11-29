--############################################### Create cross sections
xs = {}

for m=0,6 do
    xs[tostring(m)] = PhysicsTransportXSCreate()
end

-- GMesh mesh
PhysicsTransportXSSet(xs["0"],CHI_XSFILE,"materials/XS_water.cxs")
PhysicsTransportXSSet(xs["1"],CHI_XSFILE,"materials/XS_UO2.cxs")
PhysicsTransportXSSet(xs["2"],CHI_XSFILE,"materials/XS_7pMOX.cxs")
PhysicsTransportXSSet(xs["3"],CHI_XSFILE,"materials/XS_guide_tube.cxs")
PhysicsTransportXSSet(xs["4"],CHI_XSFILE,"materials/XS_4_3pMOX.cxs")
PhysicsTransportXSSet(xs["5"],CHI_XSFILE,"materials/XS_8_7pMOX.cxs")
PhysicsTransportXSSet(xs["6"],CHI_XSFILE,"materials/XS_fission_chamber.cxs")

water_xs = PhysicsTransportXSGet(xs["0"])

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
