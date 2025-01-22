function PrintPoint(pt, indent)
  print(string.rep(" ", indent) .. "1 " .. tostring(pt.x))
  print(string.rep(" ", indent) .. "2 " .. tostring(pt.y))
  print(string.rep(" ", indent) .. "3 " .. tostring(pt.z))
end

function PrintTable(t, indent)
  if not indent then
    indent = 0
  end
  strform = "%" .. tostring(indent) .. "s"

  for k, v in pairs(t) do
    if type(v) == "table" then
      print(string.rep(" ", indent) .. k .. " " .. "table")
      PrintPoint(v, indent + 2)
    else
      print(string.rep(" ", indent) .. k .. " " .. tostring(v))
    end
  end
end

print("GOLD_BEGIN")
q = squad.GaussLegendreQuadrature.Create({ N = 4, verbose = true })
print("qpoints:")
PrintTable(q.qpoints, 2)
print("weights:")
PrintTable(q.weights, 2)
print()

--
q = squad.GaussChebyshevQuadrature.Create({ N = 4, verbose = true })
print("qpoints:")
PrintTable(q.qpoints, 2)
print("weights:")
PrintTable(q.weights, 2)

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
