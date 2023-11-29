--############################################### Create cross sections
xs = {}

for m=0,1 do
    xs[tostring(m)] = chiPhysicsTransportXSCreate()
end

chiPhysicsTransportXSSet(xs["0"],CHI_XSFILE,"../transport_keigen/xs_water_g2.cxs")
chiPhysicsTransportXSSet(xs["1"],CHI_XSFILE,"../transport_keigen/xs_fuel_g2.cxs")

water_xs = chiPhysicsTransportXSGet(xs["0"])
num_groups = water_xs["num_groups"]

--############################################### Create materials
materials = {}
for m=0,1 do
    key = tostring(m)
    materials[key] = PhysicsAddMaterial("Material_"..key)
    PhysicsMaterialAddProperty(key,TRANSPORT_XSECTIONS)
    chiPhysicsMaterialSetProperty(key,TRANSPORT_XSECTIONS, EXISTING,xs[key])
end
