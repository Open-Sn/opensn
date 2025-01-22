-- read a cross section file, apply a scaling factor, and printout some components to the console

-- lua function to print entries of a table
function dump(o)
  if type(o) == "table" then
    local s = "{ "
    for k, v in pairs(o) do
      if type(k) ~= "number" then
        k = '"' .. k .. '"'
      end
      s = s .. "[" .. k .. "] = " .. dump(v) .. ","
    end
    return s .. "} "
  else
    return tostring(o)
  end
end

-- Create cross sections
my_xs = {}

my_xs["fuel"] = xs.LoadFromOpenMC(
  "../../modules/linear_boltzmann_solvers/transport_keigen/u235.h5",
  "u235",
  294.0
)

-- print to console
print("chi\n", dump(my_xs["fuel"].chi))
print("sigma total\n", dump(my_xs["fuel"].sigma_t))
chi_before = my_xs["fuel"].chi[1]
sigt_before = my_xs["fuel"].sigma_t[1]

-- scaling factor
my_xs["fuel"]:SetScalingFactor(2.0)

-- print to console again
print("\n After scaling:")
print("chi\n", dump(my_xs["fuel"].chi))
print("sigma total\n", dump(my_xs["fuel"].sigma_t))

chi_after = my_xs["fuel"]["chi"][1]
sigt_after = my_xs["fuel"]["sigma_t"][1]

log.Log(LOG_0, "chi[1] before: " .. chi_before)
log.Log(LOG_0, "chi[1] after : " .. chi_after)
log.Log(LOG_0, "sigt[1] before: " .. sigt_before)
log.Log(LOG_0, "sigt[1] after : " .. sigt_after)

-- Create a material
materials = {}
materials[0] = mat.AddMaterial("Material_fuel")
materials[0]:SetTransportXSections(my_xs["fuel"])
