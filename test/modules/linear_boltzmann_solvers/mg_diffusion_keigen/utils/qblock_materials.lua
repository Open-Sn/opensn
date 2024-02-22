--############################################### Create cross sections
xs = {}

for m=0,1 do
    xs[tostring(m)] = PhysicsTransportXSCreate()
end

PhysicsTransportXSSet(xs["0"],OPENSN_XSFILE,"../transport_keigen/xs_water_g2.xs")
PhysicsTransportXSSet(xs["1"],OPENSN_XSFILE,"../transport_keigen/xs_fuel_g2.xs")

water_xs = PhysicsTransportXSGet(xs["0"])
num_groups = water_xs["num_groups"]

--############################################### Create materials
materials = {}
for m=0,1 do
    key = tostring(m)
    materials[key] = mat.AddMaterial("Material_"..key)
    mat.AddProperty(key, TRANSPORT_XSECTIONS)
    mat.SetProperty(key, TRANSPORT_XSECTIONS, EXISTING, xs[key])
end
