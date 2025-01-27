-- Create cross sections
xss = {}

for m = 0, 6 do
  xss[tostring(m)] = xs.Create()
end

xss["0"] = xs.LoadFromOpenSn("materials/XS_water.xs")
xss["1"] = xs.LoadFromOpenSn("materials/XS_UO2.xs")
xss["2"] = xs.LoadFromOpenSn("materials/XS_7pMOX.xs")
xss["3"] = xs.LoadFromOpenSn("materials/XS_guide_tube.xs")
xss["4"] = xs.LoadFromOpenSn("materials/XS_4_3pMOX.xs")
xss["5"] = xs.LoadFromOpenSn("materials/XS_8_7pMOX.xs")
xss["6"] = xs.LoadFromOpenSn("materials/XS_fission_chamber.xs")

num_groups = xss["0"].num_groups
log.Log(LOG_0, "Num groups: " .. tostring(num_groups))

-- Create materials
materials = {}
for m = 0, 6 do
  key = tostring(m)
  materials[key] = mat.AddMaterial("Material_" .. key)
  materials[key]:SetTransportXSections(xss[key])
end
