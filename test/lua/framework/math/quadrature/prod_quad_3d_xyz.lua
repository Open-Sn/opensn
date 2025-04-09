-- Test 3D slab Gauss-Legendre Chebyshev product quadrature

n_polar = 2
n_azim = 4
pquad = aquad.CreateGLCProductQuadrature3DXYZ(n_polar, n_azim)

function rad2deg(radians)
  return radians * (180 / math.pi)
end

log.Log(LOG_0, "TEST_BEGIN")
for i = 1, n_polar * n_azim do
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
