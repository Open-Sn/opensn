// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "python/lib/py_wrappers.h"
#include "framework/event_system/event.h"
#include "framework/logging/log.h"
#include "framework/post_processors/aggregate_nodal_value_post_processor.h"
#include "framework/post_processors/cell_volume_integral_post_processor.h"
#include "framework/post_processors/post_processor.h"
#include "framework/post_processors/post_processor_printer.h"
#include "framework/post_processors/solver_info_post_processor.h"
#include "framework/utils/utils.h"
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

// Wrap post processors
void
WrapPostProcessor(py::module& post)
{
  // clang-format off
  // post processor
  auto pp = py::class_<PostProcessor, std::shared_ptr<PostProcessor>>(
    post,
    "PostProcessor",
    R"(
    Base class for all post-processors.

    Wrapper of :cpp:class:`opensn::PostProcessor`.
    )"
  );
  pp.def(
    "GetValue",
    [](PostProcessor& self)
    {
      const ParameterBlock& value = self.GetValue();
      return value.GetValue<double>();
    },
    R"(
    ???
    )"
  );
  pp.def(
    "Execute",
    [](PostProcessor& self, const std::string& event_name){
      self.Execute(Event(event_name));
    },
    R"(
    ???

    Parameters
    ----------
    event_name: str, default='ManualExecutation'
        ???
    )",
    py::arg("event_name") = "ManualExecutation"
  );

  // solver info post processor
  auto solver_info_pp = py::class_<SolverInfoPostProcessor,
                                   std::shared_ptr<SolverInfoPostProcessor>,
                                   PostProcessor>(
    post,
    "SolverInfoPostProcessor",
    R"(
    Post-processor for basic info of a Solver.

    This solver's execution does not filter whether solver events are for the relevant solver. This
    is done to avoid differing time-histories.

    (really unclear, rephrasing is needed)
    ???

    Wrapper of :cpp:class:`opensn::SolverInfoPostProcessor`.
    )"
  );
  solver_info_pp.def(
    py::init(
      [](py::kwargs& params)
      {
        return SolverInfoPostProcessor::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a solver info post processor object.

    Parameters
    ----------
    name: str
        Name of the post processor; used throughout the program.
    execute_on: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        List of events at which the post-processor executes.
    print_on: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        List of events at which the post-processor prints, ensuring these events are also set for the
        post-processor printer.
    initial_value: Dict, default={}
        Dictionary of initial value.
    print_numeric_format: {'fixed', 'floating_point', 'scientific', 'general'}, default='general'
        Numeric format for printing.
    print_precision: int, default=6
        Number of digits displayed after the decimal point.
    solvername_filter: str, default=''
        Filter to control update events to execute only on calls from the specified solver.
    solver: pyopensn.solver.Solver
        Target solver.
    info: Dict, default={}
        Dictionary requiring at minimum the parameter ``name`` to pass to the solver. Each solver
        may define additional custom parameters to retrieve various types of information.
    )"
  );

  // aggregate nodal value post processor
  auto agg_nodel_value_pp = py::class_<AggregateNodalValuePostProcessor,
                                       std::shared_ptr<AggregateNodalValuePostProcessor>,
                                       PostProcessor>(
    post,
    "AggregateNodalValuePostProcessor",
    R"(
    Aggregator for nodal values.

    Gets the max/min/avg nodal value of a field function among nodal values.

    Wrapper of :cpp:class:`opensn::AggregateNodalValuePostProcessor`.
    )"
  );
  agg_nodel_value_pp.def(
    py::init(
      [](py::kwargs& params)
      {
        return AggregateNodalValuePostProcessor::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct an aggregate nodal value post processor object.

    Parameters
    ----------
    name: str
        Name of the post processor; used throughout the program.
    execute_on: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        List of events at which the post-processor executes.
    print_on: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        List of events at which the post-processor prints, ensuring these events are also set for the
        post-processor printer.
    initial_value: Dict, default={}
        Dictionary of initial value.
    print_numeric_format: {'fixed', 'floating_point', 'scientific', 'general'}, default='general'
        Numeric format for printing.
    print_precision: int, default=6
        Number of digits displayed after the decimal point.
    solvername_filter: str, default=''
        Filter to control update events to execute only on calls from the specified solver.
    operation: {'max', 'min', 'avg'}
        The required operation to be performed.
    )"
  );

  // cell volume integral post processor
  auto cell_volume_int_pp = py::class_<CellVolumeIntegralPostProcessor,
                                       std::shared_ptr<CellVolumeIntegralPostProcessor>,
                                       PostProcessor>(
    post,
    "CellVolumeIntegralPostProcessor",
    R"(
    Compute the volumetric integral of a field-function.
    )"
  );
  cell_volume_int_pp.def(
    py::init(
      [](py::kwargs& params)
      {
        return CellVolumeIntegralPostProcessor::Create(kwargs_to_param_block(params));
      }
    ),
    R"(
    Construct a cell volume integral post processor object.

    Parameters
    ----------
    name: str
        Name of the post processor; used throughout the program.
    execute_on: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        List of events at which the post-processor executes.
    print_on: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        List of events at which the post-processor prints, ensuring these events are also set for the
        post-processor printer.
    initial_value: Dict, default={}
        Dictionary of initial value.
    print_numeric_format: {'fixed', 'floating_point', 'scientific', 'general'}, default='general'
        Numeric format for printing.
    print_precision: int, default=6
        Number of digits displayed after the decimal point.
    solvername_filter: str, default=''
        Filter to control update events to execute only on calls from the specified solver.
    compute_volume_average: bool, default=False
        Flag, when true will compute the volume average of the post-processor.
    )"
  );
  // clang-format on
}

// Wrap printer options
void
WrapPrinter(py::module& post)
{
  // clang-format off
  // print
  post.def(
    "Print",
    [](py::sequence& pp_list)
    {
      // convert to vector of pointers
      std::vector<const PostProcessor*> cpp_pp_list;
      cpp_pp_list.reserve(pp_list.size());
      for (py::handle py_pp: pp_list)
      {
        cpp_pp_list.push_back(py_pp.cast<const PostProcessor*>());
      }
      // get printer
      PostProcessorPrinter& printer = PostProcessorPrinter::GetInstance();
      const std::string output = printer.GetPrintedPostProcessors(cpp_pp_list);
      log.Log() << output;
    },
    R"(
    Print a list of post-processors.

    Parameters
    ---------
    pp_list: List[pyopensn.post.PostProcessor]
        List of post-processors to print.
    )",
    py::arg("pp_list")
  );
  // set options
  post.def(
    "SetPrinterOptions",
    [](py::kwargs& params)
    {
      // get printer
      PostProcessorPrinter& printer = PostProcessorPrinter::GetInstance();
      for (auto [key, value] : params)
      {
        // compute hash for key
        std::string key_name = key.cast<std::string>();
        std::uint32_t key_hash = hash_djb2a(key_name);
        // set based on hash
        switch (key_hash)
        {
          case "scalar_pp_table_format"_hash:
          {
            std::string option = value.cast<std::string>();
            if (option == "vertical")
            {
              printer.SetScalarPPTableFormat(ScalarPPTableFormat::VERTICAL);
            }
            else if (option == "horizontal")
            {
              printer.SetScalarPPTableFormat(ScalarPPTableFormat::HORIZONTAL);
            }
            else
            {
              throw std::invalid_argument("Unsupported format \"" + option +
                                          "\" specified for option \"scalar_pp_table_format\"");
            }
            log.Log() << "PostProcessorPrinter scalar_pp_table_format set to " << option;
            break;
          }
          case "events_on_which_to_print_postprocs"_hash:
          {
            std::vector<std::string> events;
            py::sequence py_value = value.cast<py::sequence>();
            events.reserve(py_value.size());
            for (py::handle py_event: py_value)
            {
              events.push_back(py_event.cast<std::string>());
            }
            printer.SetEventsOnWhichPrintPPs(events);
            log.Log() << "PostProcessorPrinter events_on_which_to_print_postprocs set";
            break;
          }
          case "print_scalar_time_history"_hash:
          {
            printer.SetPrintScalarTimeHistory(value.cast<bool>());
            break;
          }
          case "print_vector_time_history"_hash:
          {
            printer.SetPrintVectorTimeHistory(value.cast<bool>());
            break;
          }
          case "per_column_size_scalars"_hash:
          {
            printer.SetScalarPerColumnSize(value.cast<bool>());
            break;
          }
          case "per_column_size_vectors"_hash:
          {
            printer.SetVectorPerColumnSize(value.cast<bool>());
            break;
          }
          case "table_column_limit"_hash:
          {
            printer.SetTableColumnLimit(value.cast<size_t>());
            break;
          }
          case "time_history_limit"_hash:
          {
            printer.SetTimeHistoryLimit(value.cast<size_t>());
            break;
          }
          case "csv_filename"_hash:
          {
            printer.SetCSVFilename(value.cast<std::string>());
            break;
          }
          default:
          {
            throw std::invalid_argument("Invalid option \"" + key_name + "\"");
          }
        }
      }
    },
    R"(
    Set printer options.

    Parameters
    ----------
    scalar_pp_table_format: {'vertical', 'horizontal'}, default='vertical'
        The table format with which to print scalar.
    events_on_which_to_print_postprocs: List[str], default=['SolverInitialized', 'SolverAdvanced', 'SolverExecuted', 'ProgramExecuted']
        A list of events on which to print post-processors.
    print_scalar_time_history: bool, default=True
        Control whether a time history of scalar post-processors are printed. If false, only the
        latest version will be printed.
    print_vector_time_history: bool, default=True
        Control whether a time history of vector post-processors are printed. If false, only the
        latest version will be printed.
    per_column_size_scalars: bool, default=True
        Control the sizing of printed columns. If false, all the columns will be the same size.
    per_column_size_vectors: bool, default=True
        Control the sizing of printed columns. If false, all the columns will be the same size.
    table_column_limit: int, default=120
        The maximum column, if reached, would cause tables to be wrapped. A minimum limit of 80 is
        automatically enforced.
    time_history_limit: int, default=15
        Maximum amount of time values to show in post-processor histories. A maximum of 1000 is
        automatically enforced.
    csv_filename: str, default=''
        If not empty, a file will be printed with all the post-processors formatted as comma
        seperated values.
    )"
  );
  // clang-format on
}

// Wrap the post-processing components of OpenSn.
void
py_post(py::module& pyopensn)
{
  py::module post = pyopensn.def_submodule("post", "Post-processing module.");
  WrapPostProcessor(post);
  WrapPrinter(post);
}

} // namespace opensn
