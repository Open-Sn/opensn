/**\defgroup doc_PostProcessors H Post-Processors
*
* # 1 What is a Post-Processor?
* A post-processor, implemented via `chi::PostProcessor`, computes simple
* quantities derived from simulation entities (i.e., meshers, solvers, etc.).
* The most common form of a post-processor is a `SCALAR` which computes a single
* value, however, there is a also a `VECTOR` type and an `ARBITRARY` type.
*
* # 2 How are post-processors created
* Below is an example of the simplest (and first created) post-processors called
* the `chi::SolverInfoPostProcessor`. This post-processor obtains information
* from a solver, where each solver can override the virtual method
* `chi_physics::Solver::GetInfo` to provide any piece(s) of information. In this
* case the `prk::TransientSolver` is a mesh-less solver that computes the
* time-dependent neutron population of a Point-Reactor model.
\code
-- Example Point-Reactor Kinetics solver
phys0 = prk.TransientSolver.Create({ initial_source = 0.0 })

pp0 = chi.SolverInfoPostProcessor.Create
({
  name = "neutron_population",
  solver = phys0,
  info = {name = "neutron_population"},
  print_on = { "ProgramExecuted" }
})
\endcode
Post-processors have a few standard parameters.The `name` parameter is always
required among all post-processors. Additionally we have the following:
- `execute_on` controlling the events on which to execute the post-processor
  (more on this later).
- `print_on` controlling the events on which to print the post-processor.
- `initial_value` can be used to set an initial value for the post-processor.
- other commands to control various behaviors see \ref chi__PostProcessor.
*
* # 3 How are post-processors executed?
Post-processors are executed using the event system (see \ref doc_EventSystem).
Each post-processor can be subscribed to objects deriving from
`chi::EventPublisher` that call the Post-Processors at different stages
(Subscriber design pattern).
One example is the physics solver-system which has a wrapper to call the post
processors on calls to `Initialize`, `Execute`, `Step` and `Advance`

## 3.1 Default `execute_on` events
By default a post-processor subscribes to the following events for execution:
- `SolverInitialized`
- `SolverAdvanced`
- `SolverExecuted`
- `ProgramExecuted`

This behavior can be customized by using the `execute_on` parameter. For
example: \code pp0 = chi.SolverInfoPostProcessor.Create
({
  name = "neutron_population",
  solver = phys0,
  info = {name = "neutron_population"},
  execute_on = {"ProgramExecuted"},
  print_on = { "ProgramExecuted" }
})
\endcode

## 3.2 Post-Processor output controls
Post-Processors are printed to console or file using the
`chi::PostProcessorPrinter` singleton. It has the options described in
\ref chi__PostProcessorPrinterOptions which can be set using
\ref chi__PostProcessorPrinterSetOptions.

### 3.2.1 Printing to console output
The following example shows the printing of a simple post processor.
\code
-- Example Point-Reactor Kinetics solver
phys0 = prk.TransientSolver.Create({ initial_source = 0.0 })

pp0 = chi.SolverInfoPostProcessor.Create
({
  name = "neutron_population",
  solver = phys0,
  info = {name = "neutron_population"},
  print_on = { "ProgramExecuted" }
})

SolverInitialize(phys0)

for t=1,20 do
  SolverStep(phys0)
  time = chiPRKGetParam(phys0, "time_next")
  print(t, time,
        chiPRKGetParam(phys0, "population_next"),
        chiPRKGetParam(phys0, "period"))

  SolverAdvance(phys0)
  if (time > 0.1) then
    prk.SetParam(phys0, "rho", 0.8)
  end
end
\endcode
for which the last portion of the output would be:
\verbatim
[0]  Final program time 00:00:00
[0]  2023-09-04 08:32:02 ChiTech finished execution of scratch/RPK/rpk1.lua
[0]
[0]  SCALAR post-processors history at event "ProgramExecuted"
[0]  *----------*--------------------*
[0]  | Time     | neutron_population |
[0]  *----------*--------------------*
[0]  | 0.060000 |           1.000000 |
[0]  | 0.070000 |           1.000000 |
[0]  | 0.080000 |           1.000000 |
[0]  | 0.090000 |           1.000000 |
[0]  | 0.100000 |           1.000000 |
[0]  | 0.110000 |           1.000000 |
[0]  | 0.120000 |           3.285132 |
[0]  | 0.130000 |           4.316595 |
[0]  | 0.140000 |           4.807549 |
[0]  | 0.150000 |           5.065645 |
[0]  | 0.160000 |           5.223603 |
[0]  | 0.170000 |           5.338682 |
[0]  | 0.180000 |           5.435592 |
[0]  | 0.190000 |           5.524998 |
[0]  | 0.200000 |           5.611508 |
[0]  |   Latest |           5.611508 |
[0]  *----------*--------------------*
[0]
[0]  SCALAR post-processors latest values at event "ProgramExecuted"
[0]  *------------------------------*-----------------*
[0]  | Post-Processor Name          | Value           |
[0]  *------------------------------*-----------------*
[0]  | neutron_population(latest)   |        5.611508 |
[0]  *------------------------------*-----------------*
\endverbatim

Here we can see that the `chi::PostProcessorPrinter` printed both the time
history of the post-processor as well its latest value.

The time history can be suppressed by adding the following to the input
\code
chi.PostProcessorPrinterSetOptions
({
  print_scalar_time_history = false
})
\endcode
which results in only the latest value being printed
\verbatim
[0]  Final program time 00:00:00
[0]  2023-09-04 08:34:41 ChiTech finished execution of scratch/RPK/rpk1.lua
[0]
[0]  SCALAR post-processors latest values at event "ProgramExecuted"
[0]  *------------------------------*-----------------*
[0]  | Post-Processor Name          | Value           |
[0]  *------------------------------*-----------------*
[0]  | neutron_population(latest)   |        5.611508 |
[0]  *------------------------------*-----------------*
\endverbatim

### 3.2.2 Exporting to CSV
The time history of a post-processor can also be exported to a
Comma Separated Value file (CSV-file) by setting the `csv_filename` option in
\ref chi__PostProcessorPrinterSetOptions. Example:
\code
chi.PostProcessorPrinterSetOptions
({
  csv_filename = "rpk1.csv"
})
\endcode

### 3.2.3 Manual printing of post-processors
Post-processors can be printing manually whenever needed by using either a list
of post-processor names or handles. For example:
\code
chi.PrintPostProcessors({"test_arb", "period(s)6"}) --using names
chi.PrintPostProcessors({pp0, pp1})                 --using handles
\endcode

### 3.2.4 Controlling numeric formats
Each post-processor has the options `print_numeric_format` and `print_precision`
that controls how numbers are printed. For more about this see
\ref chi__PostProcessor.
*
* */
