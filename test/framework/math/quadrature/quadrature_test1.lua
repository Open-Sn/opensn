function PrintTable(t, indent)
  if not indent then indent = 0 end
  strform = "%"..tostring(indent).."s"

  for k,v in pairs(t) do
    if (type(v) == "table") then
      print(string.rep(" ", indent)..k.." ".."table")
      PrintTable(v, indent+2)
    else
      print(string.rep(" ", indent)..k.." "..tostring(v))
    end
  end
end

print("GOLD_BEGIN")
q = chi_math.QuadratureGaussLegendre.Create({N = 4, verbose = true})

qdata = chi_math.Get1DQuadratureData(q)

print("qpoints:")
PrintTable(qdata.qpoints, 2)
print("weights:")
PrintTable(qdata.weights, 2)
print()

--################################################
q = chi_math.QuadratureGaussChebyshev.Create({N = 4, verbose = true})

qdata = chi_math.Get1DQuadratureData(q)

print("qpoints:")
PrintTable(qdata.qpoints, 2)
print("weights:")
PrintTable(qdata.weights, 2)

print("chiLegendre(0, 0.25)", chiLegendre(0, 0.25))
print("chiLegendre(1, 0.25)", chiLegendre(1, 0.25))
print("chiLegendreDerivative(0, 0.25)", chiLegendreDerivative(0, 0.25))
print("chiLegendreDerivative(1, 0.25)", chiLegendreDerivative(1, 0.25))

print("chiYlm(0, 0, 45*math.pi/180.0, 45*math.pi/180.0)", chiYlm(0, 0, 45*math.pi/180.0, 45*math.pi/180.0))
print("chiYlm(1, 0, 45*math.pi/180.0, 45*math.pi/180.0)", chiYlm(1, 0, 45*math.pi/180.0, 45*math.pi/180.0))

print("GOLD_END")