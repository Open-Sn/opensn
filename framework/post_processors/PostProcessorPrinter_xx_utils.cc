#include "PostProcessorPrinter.h"

#include "post_processors/PostProcessor.h"
#include "event_system/Event.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <string>

namespace chi
{

// ##################################################################
std::vector<const PostProcessor*>
PostProcessorPrinter::GetScalarPostProcessorsList(const Event& event)
{
  std::vector<const PostProcessor*> scalar_pp_list;
  for (const auto& pp : Chi::postprocessor_stack)
  {
    const auto& scope = pp->PrintScope();

    // Check whether the pp wants to be printed on this event
    if (std::find(scope.begin(), scope.end(), event.Name()) == scope.end()) continue;

    if (pp->Type() == PPType::SCALAR) scalar_pp_list.push_back(&(*pp));
  }

  return scalar_pp_list;
}

// ##################################################################
std::vector<const PostProcessor*>
PostProcessorPrinter::GetVectorPostProcessorsList(const Event& event)
{
  std::vector<const PostProcessor*> scalar_pp_list;
  for (const auto& pp : Chi::postprocessor_stack)
  {
    const auto& scope = pp->PrintScope();

    // Check whether the pp wants to be printed on this event
    if (std::find(scope.begin(), scope.end(), event.Name()) == scope.end()) continue;

    if (pp->Type() == PPType::VECTOR) scalar_pp_list.push_back(&(*pp));
  }

  return scalar_pp_list;
}

// ##################################################################
std::vector<const PostProcessor*>
PostProcessorPrinter::GetArbitraryPostProcessorsList(const Event& event)
{
  std::vector<const PostProcessor*> scalar_pp_list;
  for (const auto& pp : Chi::postprocessor_stack)
  {
    const auto& scope = pp->PrintScope();

    // Check whether the pp wants to be printed on this event
    if (std::find(scope.begin(), scope.end(), event.Name()) == scope.end()) continue;

    if (pp->Type() == PPType::ARBITRARY) scalar_pp_list.push_back(&(*pp));
  }

  return scalar_pp_list;
}

// ##################################################################
std::vector<std::vector<std::string>>
PostProcessorPrinter::BuildPPHistoryMatrix(size_t timehistsize,
                                           size_t time_history_limit,
                                           const std::vector<const PostProcessor*>& pp_sub_list)
{
  if (pp_sub_list.empty()) return {};

  //+2 top header + bottom header
  const size_t num_rows = std::min(size_t(time_history_limit + 2), timehistsize + 2);
  const size_t num_cols = pp_sub_list.size() + 1; //+1 time column
  const size_t offset = std::max(0, int(timehistsize) - int(time_history_limit));

  const auto& front_time_hist = pp_sub_list.front()->GetTimeHistory();

  typedef std::vector<std::string> VecStr;
  typedef std::vector<VecStr> MatStr;
  MatStr value_matrix(num_rows, VecStr(num_cols, ""));

  // Do the header first
  value_matrix[0][0] = "Time";
  for (size_t j = 1; j <= pp_sub_list.size(); ++j)
    value_matrix[0][j] = pp_sub_list.at(j - 1)->Name();

  // Now the time values
  for (size_t t = 0; t < (num_rows - 2); ++t)
  {
    for (size_t j = 0; j <= pp_sub_list.size(); ++j)
    {
      if (j == 0) value_matrix[t + 1][j] = std::to_string(front_time_hist[t + offset].time_);
      else
      {
        const auto& pp = pp_sub_list.at(j - 1);
        value_matrix[t + 1][j] =
          pp->ConvertValueToString(pp->GetTimeHistory().at(t + offset).value_);
      }
    } // for j
  }   // for t

  // Now the last row
  {
    size_t t = num_rows - 1;
    value_matrix[t][0] = "Latest";
    for (size_t j = 0; j < pp_sub_list.size(); ++j)
    {
      const auto& pp = pp_sub_list.at(j);
      value_matrix[t][j + 1] = pp->ConvertValueToString(pp->GetValue());
    }
  }

  return value_matrix;
}

} // namespace chi