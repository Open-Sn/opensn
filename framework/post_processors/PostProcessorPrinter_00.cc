#include "PostProcessorPrinter.h"
#include "event_system/SystemWideEventPublisher.h"
#include "event_system/EventSubscriber.h"
#include "event_system/Event.h"

#include "post_processors/PostProcessor.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>

/**Small utility macro for joining two words.*/
#define JoinWordsA(x, y) x##y
/**IDK why this is needed. Seems like counter doesnt work properly without it*/
#define JoinWordsB(x, y) JoinWordsA(x, y)

std::shared_ptr<chi::PPPrinterSubscribeHelper>
  chi::PostProcessorPrinter::helper_ptr_ =
    std::make_shared<PPPrinterSubscribeHelper>(
      PostProcessorPrinter::GetInstance());

static char JoinWordsB(unique_var_name_ppp_, __COUNTER__) =
  chi::PostProcessorPrinter::SubscribeToSystemWideEventPublisher();

namespace chi
{

// ##################################################################
PPPrinterSubscribeHelper::PPPrinterSubscribeHelper(
  PostProcessorPrinter& printer_ref)
  : printer_ref_(printer_ref)
{
}

// ##################################################################
/**Overrides base class to forward the call to the printer.*/
void PPPrinterSubscribeHelper::ReceiveEventUpdate(const Event& event)
{
  printer_ref_.ReceiveEventUpdate(event);
}

// ##################################################################
PostProcessorPrinter::PostProcessorPrinter()
  : events_on_which_to_print_postprocs_({"SolverInitialized",
                                         "SolverAdvanced",
                                         "SolverExecuted",
                                         "ProgramExecuted"})
{
}

// ##################################################################
PostProcessorPrinter& PostProcessorPrinter::GetInstance()
{
  static PostProcessorPrinter instance;

  return instance;
}

// ##################################################################
char PostProcessorPrinter::SubscribeToSystemWideEventPublisher()
{
  auto helper_ptr = PostProcessorPrinter::helper_ptr_;

  auto& publisher = SystemWideEventPublisher::GetInstance();
  auto subscriber_ptr = std::dynamic_pointer_cast<EventSubscriber>(helper_ptr);

  ChiLogicalErrorIf(
    not subscriber_ptr,
    "Failure to cast chi::PPPrinterSubscribeHelper to chi::EventSubscriber");

  publisher.AddSubscriber(subscriber_ptr);

  return 0;
}

// ##################################################################
void PostProcessorPrinter::SetScalarPPTableFormat(ScalarPPTableFormat format)
{
  scalar_pp_table_format_ = format;
}

// ##################################################################
void PostProcessorPrinter::SetEventsOnWhichPrintPPs(
  const std::vector<std::string>& events)
{
  events_on_which_to_print_postprocs_ = events;
}

// ##################################################################
void PostProcessorPrinter::SetPrintScalarTimeHistory(bool value)
{
  print_scalar_time_history_ = value;
}

// ##################################################################
void PostProcessorPrinter::SetPrintVectorTimeHistory(bool value)
{
  print_vector_time_history_ = value;
}

// ##################################################################
void PostProcessorPrinter::SetScalarPerColumnSize(bool value)
{
  per_column_size_scalars_ = value;
}

// ##################################################################
void PostProcessorPrinter::SetVectorPerColumnSize(bool value)
{
  per_column_size_vectors_ = value;
}

// ##################################################################
void PostProcessorPrinter::SetTableColumnLimit(size_t limit)
{
  table_column_limit_ = std::max(limit, size_t(80));
}

// ##################################################################
void PostProcessorPrinter::SetTimeHistoryLimit(size_t limit)
{
  time_history_limit_ = std::min(limit, size_t(1000));
}

// ##################################################################
void PostProcessorPrinter::SetCSVFilename(const std::string& csv_filename)
{
  csv_filename_ = csv_filename;
}

// ##################################################################
void PostProcessorPrinter::ReceiveEventUpdate(const Event& event)
{
  {
    auto& vec = events_on_which_to_print_postprocs_;
    auto it = std::find(vec.begin(), vec.end(), event.Name());
    if (it != vec.end()) PrintPostProcessors(event);
  }
}

// ##################################################################
void PostProcessorPrinter::PrintPostProcessors(const Event& event) const
{
  const auto scalar_pps = GetScalarPostProcessorsList(event);
  {
    if (not print_scalar_time_history_)
      PrintPPsLatestValuesOnly("SCALAR", scalar_pps, event);
    else
      PrintPPsTimeHistory(
        "SCALAR", scalar_pps, event, per_column_size_scalars_);

    // If we are not printing the latest values, then how would we get values
    // suitable for regression tests. This is how.
    if (print_scalar_time_history_ and event.Name() == "ProgramExecuted")
      PrintPPsLatestValuesOnly("SCALAR", scalar_pps, event);
  }

  const auto vector_pps = GetVectorPostProcessorsList(event);
  {
    if (not print_vector_time_history_)
      PrintPPsLatestValuesOnly("VECTOR", vector_pps, event);
    else
      PrintPPsTimeHistory(
        "VECTOR", vector_pps, event, per_column_size_vectors_);

    // If we are not printing the latest values, then how would we get values
    // suitable for regression tests. This is how.
    if (print_vector_time_history_ and event.Name() == "ProgramExecuted")
      PrintPPsLatestValuesOnly("VECTOR", vector_pps, event);
  }

  if (not csv_filename_.empty() and event.Name() == "ProgramExecuted")
    PrintCSVFile(event);
}

// ##################################################################
/**A manual means to print a post processor.*/
std::string PostProcessorPrinter::GetPrintedPostProcessors(
  const std::vector<const PostProcessor*>& pp_list) const
{
  std::stringstream outstr;

  typedef std::pair<std::string, std::string> PPNameAndVal;
  std::vector<PPNameAndVal> scalar_ppnames_and_vals;
  for (const auto& pp : pp_list)
  {
    const auto& value = pp->GetValue();
    const auto value_str = pp->ConvertValueToString(value);

    scalar_ppnames_and_vals.emplace_back(pp->Name(), value_str);
  } // for pp

  if (not scalar_ppnames_and_vals.empty())
  {
    if (scalar_pp_table_format_ == ScalarPPTableFormat::HORIZONTAL)
      outstr << PrintPPsHorizontal(scalar_ppnames_and_vals, 0);
    else if (scalar_pp_table_format_ == ScalarPPTableFormat::VERTICAL)
      outstr << PrintPPsVertical(scalar_ppnames_and_vals, 0);
  }

  return outstr.str();
}

} // namespace chi