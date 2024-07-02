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

my_xs["fuel"] = xs.Create()
xs.Set(
  my_xs["fuel"],
  OPENMC_XSLIB,
  "../../modules/linear_boltzmann_solvers/transport_keigen/u235.h5",
  294.0,
  "u235"
)

-- print to console
read_xs = xs.Get(my_xs["fuel"])
print("chi\n", dump(read_xs["chi"]))
print("sigma total\n", dump(read_xs["sigma_t"]))
chi_before = read_xs["chi"][1]
sigt_before = read_xs["sigma_t"][1]

-- scaling factor
xs.SetScalingFactor(my_xs["fuel"], 2.0)

-- print to console again
read_xs = xs.Get(my_xs["fuel"])
print("\n After scaling:")
print("chi\n", dump(read_xs["chi"]))
print("sigma total\n", dump(read_xs["sigma_t"]))

chi_after = read_xs["chi"][1]
sigt_after = read_xs["sigma_t"][1]

log.Log(LOG_0, "chi[1] before: " .. chi_before)
log.Log(LOG_0, "chi[1] after : " .. chi_after)
log.Log(LOG_0, "sigt[1] before: " .. sigt_before)
log.Log(LOG_0, "sigt[1] after : " .. sigt_after)

-- Create a material
materials = {}
materials[0] = mat.AddMaterial("Material_fuel")
mat.SetProperty(materials[0], TRANSPORT_XSECTIONS, EXISTING, my_xs["fuel"])
