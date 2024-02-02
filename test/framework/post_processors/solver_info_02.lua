-- Post-Processor test with lots of post-processors
-- Testing table wrapping and getting the value of a post-processor by both
-- handle and name

-- Example Point-Reactor Kinetics solver
phys0 = prk.TransientSolver.Create({ initial_source = 0.0 })

for k = 1, 20 do
  SolverInfoPostProcessor.Create({
    name = "neutron_population" .. tostring(k),
    solver = phys0,
    info = { name = "neutron_population" },
    print_on = { "" }
  })
end
pp21 = SolverInfoPostProcessor.Create({
  name = "neutron_population" .. tostring(21),
  solver = phys0,
  info = { name = "neutron_population" },
  print_on = { "" }
})

PostProcessorPrinterSetOptions({
  time_history_limit = 5,
})

SolverInitialize(phys0)

for t = 1, 20 do
  SolverStep(phys0)
  time = SolverGetInfo(phys0, "time_next")
  print(t, string.format("%.3f %.5f",time, SolverGetInfo(phys0, "population_next")))

  SolverAdvance(phys0)
  if (time > 0.1) then
    prk.SetParam(phys0, "rho", 0.8)
  end
end

print("Manual neutron_population1=",
  string.format("%.5f", PostProcessorGetValue("neutron_population1")))
print("Manual neutron_population1=",
  string.format("%.5f", PostProcessorGetValue(pp21)))
