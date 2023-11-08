#include "framework/post_processors/PostProcessorPrinter.h"

#include "framework/parameters/input_parameters.h"
#include "framework/utils/chi_utils.h"

#include "framework/console/chi_console.h"
#include "framework/ChiObjectFactory.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace chi
{

InputParameters PostProcessorPrinterOptions();
InputParameters GetSyntax_PPPrinterSetOptions();
ParameterBlock PostProcessorPrinterSetOptions(const InputParameters& params);

RegisterSyntaxBlock(chi, PostProcessorPrinterOptions, PostProcessorPrinterOptions);

RegisterWrapperFunction(/*namespace_in_lua=*/chi,
                        /*name_in_lua=*/PostProcessorPrinterSetOptions,
                        /*syntax_function=*/GetSyntax_PPPrinterSetOptions,
                        /*actual_function=*/PostProcessorPrinterSetOptions);

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

  using namespace chi_data_types;
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

InputParameters
GetSyntax_PPPrinterSetOptions()
{
  InputParameters params;

  params.SetGeneralDescription(
    "Wrapper function to set options of the singleton PostProcessorPrinter.");
  params.SetDocGroup("doc_PPUtils");

  params.AddRequiredParameterBlock("arg0", "Options parameter block");
  params.LinkParameterToBlock("arg0", "chi::PostProcessorPrinterOptions");

  return params;
}

ParameterBlock
PostProcessorPrinterSetOptions(const InputParameters& params)
{
  auto& printer = PostProcessorPrinter::GetInstance();

  const auto& set_params = params.ParametersAtAssignment().GetParam("arg0");

  for (const auto& param : set_params)
  {
    const std::string param_name = param.Name();

    uint32_t param_name_hash = chi::hash_djb2a(param_name);
    // clang-format off
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
          ChiInvalidArgument(
            "Unsupported format \"" + option +
            "\" specified for option \"scalar_pp_table_format\"");

        Chi::log.Log() << "PostProcessorPrinter scalar_pp_table_format set to "
                       << option;
        break;
      }
      case "events_on_which_to_print_postprocs"_hash:
      {
        const auto list = param.GetVectorValue<std::string>();

        printer.SetEventsOnWhichPrintPPs(list);

        Chi::log.Log()
          << "PostProcessorPrinter events_on_which_to_print_postprocs set";
        break;
      }
      case "print_scalar_time_history"_hash:
        { printer.SetPrintScalarTimeHistory(param.GetValue<bool>()); break; }
      case "print_vector_time_history"_hash:
        { printer.SetPrintVectorTimeHistory(param.GetValue<bool>()); break;}
      case "per_column_size_scalars"_hash:
        { printer.SetScalarPerColumnSize(param.GetValue<bool>()); break;}
      case "per_column_size_vectors"_hash:
        { printer.SetVectorPerColumnSize(param.GetValue<bool>()); break;}
      case "table_column_limit"_hash:
        { printer.SetTableColumnLimit(param.GetValue<size_t>()); break;}
      case "time_history_limit"_hash:
        { printer.SetTimeHistoryLimit(param.GetValue<size_t>()); break;}
      case "csv_filename"_hash:
        { printer.SetCSVFilename(param.GetValue<std::string>()); break;}
      default: ChiInvalidArgument("Invalid option \"" + param.Name() + "\"");
    }// switch
    // clang-format on
  }
  return ParameterBlock{};
}

} // namespace chi
