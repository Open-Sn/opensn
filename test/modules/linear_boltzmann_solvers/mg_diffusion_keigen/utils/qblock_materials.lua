-- Create cross sections
xss = {}

xss["0"] = xs.LoadFromOpenSn("../transport_keigen/xs_water_g2.xs")
xss["1"] = xs.LoadFromOpenSn("../transport_keigen/xs_fuel_g2.xs")

num_groups = xss["0"].num_groups

-- Create materials
materials = {}
for m = 0, 1 do
  key = tostring(m)
  materials[key] = mat.AddMaterial("Material_" .. key)
  materials[key]:SetTransportXSections(xss[key])
end
