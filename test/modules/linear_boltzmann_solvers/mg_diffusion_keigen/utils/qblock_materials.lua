--############################################### Create cross sections
xss = {}

for m = 0, 1 do
  xss[tostring(m)] = xs.Create()
end

xs.Set(xss["0"], OPENSN_XSFILE, "../transport_keigen/xs_water_g2.xs")
xs.Set(xss["1"], OPENSN_XSFILE, "../transport_keigen/xs_fuel_g2.xs")

water_xs = xs.Get(xss["0"])
num_groups = water_xs["num_groups"]

--############################################### Create materials
materials = {}
for m = 0, 1 do
  key = tostring(m)
  materials[key] = mat.AddMaterial("Material_" .. key)
  mat.SetProperty(materials[key], TRANSPORT_XSECTIONS, EXISTING, xss[key])
end
