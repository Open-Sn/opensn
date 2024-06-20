--############################################### Create cross sections
xss = {}

for m = 0, 6 do
  xss[tostring(m)] = xs.Create()
end

-- GMesh mesh
xs.Set(xss["0"], OPENSN_XSFILE, "materials/XS_water.xs")
xs.Set(xss["1"], OPENSN_XSFILE, "materials/XS_UO2.xs")
xs.Set(xss["2"], OPENSN_XSFILE, "materials/XS_7pMOX.xs")
xs.Set(xss["3"], OPENSN_XSFILE, "materials/XS_guide_tube.xs")
xs.Set(xss["4"], OPENSN_XSFILE, "materials/XS_4_3pMOX.xs")
xs.Set(xss["5"], OPENSN_XSFILE, "materials/XS_8_7pMOX.xs")
xs.Set(xss["6"], OPENSN_XSFILE, "materials/XS_fission_chamber.xs")

water_xs = xs.Get(xss["0"])

num_groups = water_xs["num_groups"]
log.Log(LOG_0, "Num groups: " .. tostring(num_groups))

--############################################### Create materials
materials = {}
for m = 0, 6 do
  key = tostring(m)
  materials[key] = mat.AddMaterial("Material_" .. key)
  mat.SetProperty(materials[key], TRANSPORT_XSECTIONS, EXISTING, xss[key])
end
