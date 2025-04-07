-- Test 1D slab Gauss-Legendre product quadrature

n_polar = 4
pquad = aquad.CreateGLProductQuadrature1DSlab(n_polar)

function rad2deg(radians)
  return radians * (180 / math.pi)
end

log.Log(LOG_0, "TEST_BEGIN")
for i = 1, n_polar do
  log.Log(
    LOG_0,
    string.format(
      "%5d | %6.4f | %7.4f | %7.3f",
      i,
      pquad.weights[i],
      rad2deg(pquad.abscissae[i].theta),
      rad2deg(pquad.abscissae[i].phi)
    )
  )
end
log.Log(LOG_0, "TEST_END")
