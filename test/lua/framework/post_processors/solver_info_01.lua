-- Post-Processor test with basic information
-- Also tests the .csv output and manual printing

-- Example Point-Reactor Kinetics solver
phys0 = prk.PRKSolver.Create({ initial_source = 0.0 })

pp0 = post.SolverInfoPostProcessor.Create({
  name = "neutron_population",
  solver = phys0,
  info = { name = "neutron_population" },
  print_on = { "ProgramExecuted" },
})

post.SetPrinterOptions({
  csv_filename = "solver_info_01.csv",
})

solver.Initialize(phys0)

for t = 1, 20 do
  solver.Step(phys0)
  time = phys0:TimeNew()
  print(t, string.format("%.3f %.5f", time, phys0:PopulationNew()))

  solver.Advance(phys0)
  if time > 0.1 then
    phys0:SetRho(0.8)
  end
end

print("Manually printing Post-Processor:")
post.Print({ pp0 })
