#include "framework/post_processors/post_processor.h"

#include "framework/console/console.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

InputParameters GetSyntax_PostProcessorGetValue();
ParameterBlock PostProcessorGetValue(const InputParameters& params);

RegisterWrapperFunction(chi,
                        PostProcessorGetValue,
                        GetSyntax_PostProcessorGetValue,
                        PostProcessorGetValue);

InputParameters
GetSyntax_PostProcessorGetValue()
{
  InputParameters params;

  params.SetGeneralDescription(
    "Wrapper function to retrieve the current value of a post-processor");
  params.SetDocGroup("doc_PPUtils");

  params.AddRequiredParameter<size_t>("arg0", "Handle or name of the post-processor");
  params.SetParameterTypeMismatchAllowed("arg0");

  return params;
}

ParameterBlock
PostProcessorGetValue(const InputParameters& params)
{
  const auto& param = params.GetParam("arg0");
  if (param.Type() == ParameterBlockType::STRING)
  {
    const auto pp_name = param.GetValue<std::string>();

    for (const auto& pp_ptr : Chi::postprocessor_stack)
      if (pp_ptr->Name() == pp_name) return pp_ptr->GetValue();

    // If we haven't returned here
    ChiInvalidArgument("Post-processor with name \"" + pp_name + "\" not found.");
  }
  else if (param.Type() == ParameterBlockType::INTEGER)
  {
    const auto pp_handle = param.GetValue<size_t>();
    const auto& pp =
      Chi::GetStackItem<PostProcessor>(Chi::postprocessor_stack, pp_handle, __FUNCTION__);

    return pp.GetValue();
  }
  else
    ChiInvalidArgument("Accepts only STRING or INTEGER for arg0.");

  return ParameterBlock{};
}
} // namespace opensn
