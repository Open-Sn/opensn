#include "framework/post_processors/PostProcessorPrinter.h"
#include "framework/post_processors/PostProcessor.h"

#include "framework/console/chi_console.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace chi
{

InputParameters GetSyntax_PrintPostProcessors();
ParameterBlock PrintPostProcessors(const InputParameters& params);

RegisterWrapperFunction(/*namespace_in_lua=*/chi,
                        /*name_in_lua=*/PrintPostProcessors,
                        /*syntax_function=*/GetSyntax_PrintPostProcessors,
                        /*actual_function=*/PrintPostProcessors);

InputParameters
GetSyntax_PrintPostProcessors()
{
  InputParameters params;

  params.SetGeneralDescription("Wrapper function to manually print post-processors via the "
                               "PostProcessorPrinter singleton.");
  params.SetDocGroup("doc_PPUtils");

  params.AddRequiredParameterArray("arg0", "A list of post-processor names or handles.");

  return params;
}

ParameterBlock
PrintPostProcessors(const InputParameters& params)
{
  const auto& arg_array = params.GetParam("arg0");

  arg_array.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  ChiInvalidArgumentIf(arg_array.NumParameters() == 0, "Empty array passed.");
  const auto& first_param = arg_array.GetParam(0);
  const auto first_param_type = first_param.Type();

  std::vector<const PostProcessor*> pp_list;

  //=================================== List of names supplied
  if (first_param_type == ParameterBlockType::STRING)
  {
    const auto name_list = arg_array.GetVectorValue<std::string>();

    for (const auto& name : name_list)
    {
      bool found = false;

      for (const auto& pp_smart_ptr : Chi::postprocessor_stack)
        if (pp_smart_ptr->Name() == name)
        {
          found = true;
          pp_list.push_back(&(*pp_smart_ptr));
          break;
        }

      ChiInvalidArgumentIf(not found,
                           "Post processor with name \"" + name +
                             "\" not found in the stack of post-processors");
    }
  }
  //=================================== List of handles supplied
  else if (first_param_type == ParameterBlockType::INTEGER)
  {
    const auto handle_list = arg_array.GetVectorValue<size_t>();

    for (const size_t handle : handle_list)
    {
      auto& pp = Chi::GetStackItem<PostProcessor>(Chi::postprocessor_stack, handle, __FUNCTION__);
      pp_list.push_back(&pp);
    }
  }
  else
    ChiInvalidArgument("The array is of type ARRAY<" + ParameterBlockTypeName(first_param_type) +
                       ">. Only ARRAY<STRING> or ARRAY<INTEGER> is allowed.");

  auto& printer = PostProcessorPrinter::GetInstance();

  const std::string output = printer.GetPrintedPostProcessors(pp_list);
  Chi::log.Log() << output;

  return ParameterBlock{};
}

} // namespace chi
