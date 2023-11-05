#include "opensn/framework/post_processors/PostProcessor.h"

#include "opensn/framework/console/chi_console.h"

#include "opensn/framework/event_system/Event.h"

namespace chi
{

InputParameters GetSyntax_ExecutePostProcessors();
ParameterBlock ExecutePostProcessors(const InputParameters& params);

RegisterWrapperFunction(/*namespace_in_lua=*/chi,
                        /*name_in_lua=*/ExecutePostProcessors,
                        /*syntax_function=*/GetSyntax_ExecutePostProcessors,
                        /*actual_function=*/ExecutePostProcessors);

InputParameters
GetSyntax_ExecutePostProcessors()
{
  InputParameters params;

  params.SetGeneralDescription("Wrapper function to manually execute post-processors");
  params.SetDocGroup("doc_PPUtils");

  params.AddRequiredParameterArray("arg0", "A list of post-processor names or handles.");

  return params;
}

ParameterBlock
ExecutePostProcessors(const InputParameters& params)
{
  const auto& arg_array = params.GetParam("arg0");

  arg_array.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  ChiInvalidArgumentIf(arg_array.NumParameters() == 0, "Empty array passed.");
  const auto& first_param = arg_array.GetParam(0);
  const auto first_param_type = first_param.Type();

  std::vector<PostProcessor*> pp_list;

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

  Event blank_event("ManualExecutation", 0);
  for (auto& pp : pp_list)
    pp->Execute(blank_event);

  return ParameterBlock{};
}

} // namespace chi
