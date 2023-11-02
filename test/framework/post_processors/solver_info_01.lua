-- Post-Processor test with basic information
-- Also tests the .csv output and manual printing

-- Example Point-Reactor Kinetics solver
phys0 = prk.TransientSolver.Create({ initial_source = 0.0 })

pp0 = chi.SolverInfoPostProcessor.Create
({
  name = "neutron_population",
  solver = phys0,
  info = {name = "neutron_population"},
  print_on = { "ProgramExecuted" }
})

chi.PostProcessorPrinterSetOptions
({
  csv_filename = "solver_info_01.csv"
})

chiSolverInitialize(phys0)

for t=1,20 do
  chiSolverStep(phys0)
  time = chiSolverGetInfo(phys0, "time_next")
  print(t, string.format("%.3f %.5f",time, chiSolverGetInfo(phys0, "population_next")))

  chiSolverAdvance(phys0)
  if (time > 0.1) then
    prk.SetParam(phys0, "rho", 0.8)
  end
end

print("Manually printing Post-Processor:")
chi.PrintPostProcessors({pp0, pp1})