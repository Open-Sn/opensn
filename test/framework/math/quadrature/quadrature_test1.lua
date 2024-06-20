function PrintTable(t, indent)
  if not indent then
    indent = 0
  end
  strform = "%" .. tostring(indent) .. "s"

  for k, v in pairs(t) do
    if type(v) == "table" then
      print(string.rep(" ", indent) .. k .. " " .. "table")
      PrintTable(v, indent + 2)
    else
      print(string.rep(" ", indent) .. k .. " " .. tostring(v))
    end
  end
end

print("GOLD_BEGIN")
q = squad.GaussLegendreQuadrature.Create({ N = 4, verbose = true })

qdata = math.Get1DQuadratureData(q)

print("qpoints:")
PrintTable(qdata.qpoints, 2)
print("weights:")
PrintTable(qdata.weights, 2)
print()

--################################################
q = squad.GaussChebyshevQuadrature.Create({ N = 4, verbose = true })

qdata = math.Get1DQuadratureData(q)

print("qpoints:")
PrintTable(qdata.qpoints, 2)
print("weights:")
PrintTable(qdata.weights, 2)

print("Legendre(0, 0.25)", aquad.Legendre(0, 0.25))
print("Legendre(1, 0.25)", aquad.Legendre(1, 0.25))
print("LegendreDerivative(0, 0.25)", aquad.LegendreDerivative(0, 0.25))
print("LegendreDerivative(1, 0.25)", aquad.LegendreDerivative(1, 0.25))

print(
  "Ylm(0, 0, 45*math.pi/180.0, 45*math.pi/180.0)",
  aquad.Ylm(0, 0, 45 * math.pi / 180.0, 45 * math.pi / 180.0)
)
print(
  "Ylm(1, 0, 45*math.pi/180.0, 45*math.pi/180.0)",
  aquad.Ylm(1, 0, 45 * math.pi / 180.0, 45 * math.pi / 180.0)
)

print("GOLD_END")
