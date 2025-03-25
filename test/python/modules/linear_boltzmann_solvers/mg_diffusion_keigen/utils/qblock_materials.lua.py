# Create cross sections
xss = {}

xss["0"] = xs.LoadFromOpenSn("+/transport_keigen/xs_water_g2.xs")
xss["1"] = xs.LoadFromOpenSn("+/transport_keigen/xs_fuel_g2.xs")

num_groups = xss["0"].num_groups
