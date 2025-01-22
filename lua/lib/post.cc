// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/console.h"
#include "framework/post_processors/post_processor.h"
#include "framework/post_processors/post_processor_printer.h"
#include "framework/parameters/input_parameters.h"
#include "framework/utils/utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/event_system/event.h"

using namespace opensn;

namespace opensnlua
{

InputParameters
PostProcessorPrinterOptions()
{
  InputParameters params;

  params.SetGeneralDescription("Options allowable for the PostProcessorPrinter");
  params.SetDocGroup("doc_PPUtils");

  params.AddOptionalParameter("scalar_pp_table_format",
                              "vertical",
                              "The table format with which to print scalar "
                              "PostProcessors");

  params.ConstrainParameterRange("scalar_pp_table_format",
                                 AllowableRangeList::New({"horizontal", "vertical"}));

  params.AddOptionalParameterArray(
    "events_on_which_to_print_postprocs",
    std::vector<std::string>{
      "SolverInitialized", "SolverAdvanced", "SolverExecuted", "ProgramExecuted"},
    "A list of events on which to print post-processors");

  params.AddOptionalParameter(
    "print_scalar_time_history",
    true,
    "Controls whether a time history of scalar post-processors are printed. If "
    "false, only the latest version will be printed.");

  params.AddOptionalParameter(
    "print_vector_time_history",
    true,
    "Controls whether a time history of vector post-processors are printed. If "
    "false, only the latest version will be printed.");

  params.AddOptionalParameter("per_column_size_scalars",
                              true,
                              "Controls the sizing of printed columns. If "
                              "false all the columns will be the same size.");

  params.AddOptionalParameter("per_column_size_vectors",
                              true,
                              "Controls the sizing of printed columns. If "
                              "false all the columns will be the same size.");

  params.AddOptionalParameter("table_column_limit",
                              120,
                              "The maximum column, if reached, would cause tables to be wrapped. A "
                              "minimum limit of 80 is automatically enforced.");

  params.AddOptionalParameter(
    "time_history_limit",
    15,
    "Maximum amount of time values to show in post-processor histories. A "
    "maximum of 1000 is automatically enforced.");

  params.AddOptionalParameter("csv_filename",
                              "",
                              "If not empty, a file will be printed with all the post-processors "
                              "formatted as comma seperated values.");

  return params;
}

void
PostProcessorPrinterSetOptions(const InputParameters& params)
{
  auto& printer = PostProcessorPrinter::GetInstance();

  for (const auto& param : params)
  {
    const std::string param_name = param.Name();

    uint32_t param_name_hash = opensn::hash_djb2a(param_name);
    switch (param_name_hash)
    {
      case "scalar_pp_table_format"_hash:
      {
        const auto option = param.GetValue<std::string>();
        if (option == "vertical")
          printer.SetScalarPPTableFormat(ScalarPPTableFormat::VERTICAL);
        else if (option == "horizontal")
          printer.SetScalarPPTableFormat(ScalarPPTableFormat::HORIZONTAL);
        else
          OpenSnInvalidArgument("Unsupported format \"" + option +
                                "\" specified for option \"scalar_pp_table_format\"");

        opensn::log.Log() << "PostProcessorPrinter scalar_pp_table_format set to " << option;
        break;
      }
      case "events_on_which_to_print_postprocs"_hash:
      {
        const auto list = param.GetVectorValue<std::string>();

        printer.SetEventsOnWhichPrintPPs(list);

        opensn::log.Log() << "PostProcessorPrinter events_on_which_to_print_postprocs set";
        break;
      }
      case "print_scalar_time_history"_hash:
        printer.SetPrintScalarTimeHistory(param.GetValue<bool>());
        break;
      case "print_vector_time_history"_hash:
        printer.SetPrintVectorTimeHistory(param.GetValue<bool>());
        break;
      case "per_column_size_scalars"_hash:
        printer.SetScalarPerColumnSize(param.GetValue<bool>());
        break;
      case "per_column_size_vectors"_hash:
        printer.SetVectorPerColumnSize(param.GetValue<bool>());
        break;
      case "table_column_limit"_hash:
        printer.SetTableColumnLimit(param.GetValue<size_t>());
        break;
      case "time_history_limit"_hash:
        printer.SetTimeHistoryLimit(param.GetValue<size_t>());
        break;
      case "csv_filename"_hash:
        printer.SetCSVFilename(param.GetValue<std::string>());
        break;
      default:
        OpenSnInvalidArgument("Invalid option \"" + param.Name() + "\"");
    } // switch
  }
}

void
PrintPostProcessors(const std::vector<std::shared_ptr<PostProcessor>>& postprocessor_ptr_list)
{
  std::vector<const PostProcessor*> pp_list;
  for (auto& ptr : postprocessor_ptr_list)
    pp_list.push_back(ptr.get());

  auto& printer = PostProcessorPrinter::GetInstance();

  const std::string output = printer.GetPrintedPostProcessors(pp_list);
  opensn::log.Log() << output;
}

void
ExecutePostProcessors(
  const std::vector<std::shared_ptr<opensn::PostProcessor>>& postprocessor_ptr_list)
{
  std::vector<PostProcessor*> pp_list;
  for (auto& ptr : postprocessor_ptr_list)
    pp_list.push_back(ptr.get());

  Event blank_event("ManualExecutation");
  for (auto& pp : pp_list)
    pp->Execute(blank_event);
}

} // namespace opensnlua
