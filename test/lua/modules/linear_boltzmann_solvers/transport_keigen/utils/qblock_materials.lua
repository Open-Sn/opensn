-- Create cross sections
xss = {}
xss["0"] = xs.LoadFromOpenSn("xs_water_g2.xs")
xss["1"] = xs.LoadFromOpenSn("xs_fuel_g2.xs")

num_groups = xss["0"].num_groups

-- Create materials
xs_map = {}
for m = 0, 1 do
  key = tostring(m)
  xs_map[m + 1] = { block_ids = { m }, xs = xss[key] }
end
