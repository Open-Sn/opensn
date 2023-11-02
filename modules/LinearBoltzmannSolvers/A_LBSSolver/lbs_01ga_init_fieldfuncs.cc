#include "lbs_solver.h"

#include "physics/FieldFunction/fieldfunction_gridbased.h"

#include "chi_runtime.h"

namespace lbs
{

void
LBSSolver::InitializeFieldFunctions()
{
  using namespace chi_math;

  if (not field_functions_.empty()) return;

  //============================================= Initialize Field Functions
  //                                              for flux moments
  phi_field_functions_local_map_.clear();

  for (size_t g = 0; g < groups_.size(); ++g)
  {
    for (size_t m = 0; m < num_moments_; m++)
    {
      std::string prefix;
      if (options_.field_function_prefix_option == "prefix")
      {
        prefix = options_.field_function_prefix;
        if (not prefix.empty()) prefix += "_";
      }
      if (options_.field_function_prefix_option == "solver_name") prefix = TextName() + "_";

      char buff[100];
      snprintf(
        buff, 99, "%sphi_g%03d_m%02d", prefix.c_str(), static_cast<int>(g), static_cast<int>(m));
      const std::string text_name = std::string(buff);

      auto group_ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
        text_name,                     // Field name
        discretization_,               // Spatial discretization
        Unknown(UnknownType::SCALAR)); // Unknown/Variable

      Chi::field_function_stack.push_back(group_ff);
      field_functions_.push_back(group_ff);

      phi_field_functions_local_map_[{g, m}] = field_functions_.size() - 1;
    } // for m
  }   // for g

  //============================================= Initialize power generation
  //                                              field function
  if (options_.power_field_function_on)
  {
    std::string prefix;
    if (options_.field_function_prefix_option == "prefix")
    {
      prefix = options_.field_function_prefix;
      if (not prefix.empty()) prefix += "_";
    }
    if (options_.field_function_prefix_option == "solver_name") prefix = TextName() + "_";

    auto power_ff = std::make_shared<chi_physics::FieldFunctionGridBased>(
      prefix + "power_generation",   // Field name
      discretization_,               // Spatial discretization
      Unknown(UnknownType::SCALAR)); // Unknown/Variable

    Chi::field_function_stack.push_back(power_ff);
    field_functions_.push_back(power_ff);

    power_gen_fieldfunc_local_handle_ = field_functions_.size() - 1;
  }
}

} // namespace lbs