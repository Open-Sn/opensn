#include "framework/post_processors/post_processor_printer.h"
#include "framework/event_system/system_wide_event_publisher.h"
#include "framework/event_system/event_subscriber.h"
#include "framework/event_system/event.h"
#include "framework/post_processors/post_processor.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/utils.h"
#include <set>
#include <algorithm>

/**Small utility macro for joining two words.*/
#define JoinWordsA(x, y) x##y
/**IDK why this is needed. Seems like counter doesnt work properly without it*/
#define JoinWordsB(x, y) JoinWordsA(x, y)

namespace opensn
{

std::shared_ptr<PPPrinterSubscribeHelper> PostProcessorPrinter::helper_ptr_ =
  std::make_shared<PPPrinterSubscribeHelper>(PostProcessorPrinter::GetInstance());

static char JoinWordsB(unique_var_name_ppp_,
                       __COUNTER__) = PostProcessorPrinter::SubscribeToSystemWideEventPublisher();

PPPrinterSubscribeHelper::PPPrinterSubscribeHelper(PostProcessorPrinter& printer_ref)
  : printer_ref_(printer_ref)
{
}

void
PPPrinterSubscribeHelper::ReceiveEventUpdate(const Event& event)
{
  printer_ref_.ReceiveEventUpdate(event);
}

PostProcessorPrinter::PostProcessorPrinter()
  : events_on_which_to_print_postprocs_(
      {"SolverInitialized", "SolverAdvanced", "SolverExecuted", "ProgramExecuted"})
{
}

PostProcessorPrinter&
PostProcessorPrinter::GetInstance()
{
  static PostProcessorPrinter instance;

  return instance;
}

char
PostProcessorPrinter::SubscribeToSystemWideEventPublisher()
{
  auto helper_ptr = PostProcessorPrinter::helper_ptr_;

  auto& publisher = SystemWideEventPublisher::GetInstance();
  auto subscriber_ptr = std::dynamic_pointer_cast<EventSubscriber>(helper_ptr);

  ChiLogicalErrorIf(not subscriber_ptr,
                    "Failure to cast chi::PPPrinterSubscribeHelper to chi::EventSubscriber");

  publisher.AddSubscriber(subscriber_ptr);

  return 0;
}

void
PostProcessorPrinter::SetScalarPPTableFormat(ScalarPPTableFormat format)
{
  scalar_pp_table_format_ = format;
}

void
PostProcessorPrinter::SetEventsOnWhichPrintPPs(const std::vector<std::string>& events)
{
  events_on_which_to_print_postprocs_ = events;
}

void
PostProcessorPrinter::SetPrintScalarTimeHistory(bool value)
{
  print_scalar_time_history_ = value;
}

void
PostProcessorPrinter::SetPrintVectorTimeHistory(bool value)
{
  print_vector_time_history_ = value;
}

void
PostProcessorPrinter::SetScalarPerColumnSize(bool value)
{
  per_column_size_scalars_ = value;
}

void
PostProcessorPrinter::SetVectorPerColumnSize(bool value)
{
  per_column_size_vectors_ = value;
}

void
PostProcessorPrinter::SetTableColumnLimit(size_t limit)
{
  table_column_limit_ = std::max(limit, size_t(80));
}

void
PostProcessorPrinter::SetTimeHistoryLimit(size_t limit)
{
  time_history_limit_ = std::min(limit, size_t(1000));
}

void
PostProcessorPrinter::SetCSVFilename(const std::string& csv_filename)
{
  csv_filename_ = csv_filename;
}

void
PostProcessorPrinter::ReceiveEventUpdate(const Event& event)
{
  {
    auto& vec = events_on_which_to_print_postprocs_;
    auto it = std::find(vec.begin(), vec.end(), event.Name());
    if (it != vec.end()) PrintPostProcessors(event);
  }
}

void
PostProcessorPrinter::PrintPostProcessors(const Event& event) const
{
  const auto scalar_pps = GetScalarPostProcessorsList(event);
  {
    if (not print_scalar_time_history_) PrintPPsLatestValuesOnly("SCALAR", scalar_pps, event);
    else
      PrintPPsTimeHistory("SCALAR", scalar_pps, event, per_column_size_scalars_);

    // If we are not printing the latest values, then how would we get values
    // suitable for regression tests. This is how.
    if (print_scalar_time_history_ and event.Name() == "ProgramExecuted")
      PrintPPsLatestValuesOnly("SCALAR", scalar_pps, event);
  }

  const auto vector_pps = GetVectorPostProcessorsList(event);
  {
    if (not print_vector_time_history_) PrintPPsLatestValuesOnly("VECTOR", vector_pps, event);
    else
      PrintPPsTimeHistory("VECTOR", vector_pps, event, per_column_size_vectors_);

    // If we are not printing the latest values, then how would we get values
    // suitable for regression tests. This is how.
    if (print_vector_time_history_ and event.Name() == "ProgramExecuted")
      PrintPPsLatestValuesOnly("VECTOR", vector_pps, event);
  }

  if (not csv_filename_.empty() and event.Name() == "ProgramExecuted") PrintCSVFile(event);
}

std::string
PostProcessorPrinter::GetPrintedPostProcessors(
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

void
PostProcessorPrinter::PrintPPsLatestValuesOnly(const std::string& pps_typename,
                                               const std::vector<const PostProcessor*>& pp_list,
                                               const Event& event) const
{
  if (pp_list.empty()) return;
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
      outstr << PrintPPsHorizontal(scalar_ppnames_and_vals, event.Code());
    else if (scalar_pp_table_format_ == ScalarPPTableFormat::VERTICAL)
      outstr << PrintPPsVertical(scalar_ppnames_and_vals, event.Code());
    log.Log() << "\n"
              << pps_typename << " post-processors latest values at event \"" << event.Name()
              << "\"\n"
              << outstr.str() << "\n";
  }
}

std::string
PostProcessorPrinter::PrintPPsHorizontal(
  const std::vector<std::pair<std::string, std::string>>& scalar_ppnames_and_vals, int)
{
  std::stringstream outstr;
  const size_t num_pps = scalar_ppnames_and_vals.size();

  std::vector<size_t> col_sizes;
  col_sizes.reserve(scalar_ppnames_and_vals.size());
  for (const auto& [name, valstring] : scalar_ppnames_and_vals)
    col_sizes.push_back(std::max(name.size(), valstring.size()));

  std::stringstream header1, header2, header3, body, footer;
  for (size_t p = 0; p < num_pps; ++p)
  {
    const size_t col_size = std::max(col_sizes[p], size_t(15));
    const auto& [ppname, ppval] = scalar_ppnames_and_vals[p];
    const auto ppname2 = ppname + std::string(col_size - ppname.size(), ' ');
    const auto ppval2 = std::string(col_size - ppval.size(), ' ') + ppval;
    for (size_t c = 0; c < (col_size + 3); ++c)
    {
      if (c == 0)
      {
        header1 << "*";
        header2 << "|";
        header3 << "*";
        body << "|";
        footer << "*";
      }
      else if (c == 2)
      {
        header1 << "_";
        header2 << ppname2;
        header3 << "-";
        body << ppval2;
        footer << "-";
      }
      else if (c < (col_size + 2) and c != 1)
      {
        header1 << "_";
        header3 << "-";
        footer << "-";
      }
      else
      {
        header1 << "_";
        header2 << " ";
        header3 << "-";
        body << " ";
        footer << "-";
      }
    }
  }
  header1 << "*";
  header2 << "|";
  header3 << "*";
  body << "|";
  footer << "*";

  outstr << header1.str() << "\n";
  outstr << header2.str() << "\n";
  outstr << header3.str() << "\n";
  outstr << body.str() << "\n";
  outstr << footer.str() << "\n";

  return outstr.str();
}

std::string
PostProcessorPrinter::PrintPPsVertical(
  const std::vector<std::pair<std::string, std::string>>& scalar_ppnames_and_vals, int event_code)
{
  std::stringstream outstr;

  const size_t num_pps = scalar_ppnames_and_vals.size();

  size_t max_colsize_name = scalar_ppnames_and_vals.front().first.size();
  size_t max_colsize_val = scalar_ppnames_and_vals.front().second.size();
  for (const auto& [name, valstring] : scalar_ppnames_and_vals)
  {
    max_colsize_name = std::max(max_colsize_name, name.size() + 8);
    max_colsize_val = std::max(max_colsize_val, valstring.size());
  }
  constexpr size_t min_col_size = 15;
  max_colsize_name = std::max(max_colsize_name, min_col_size + 5);
  max_colsize_val = std::max(max_colsize_val, min_col_size);

  const std::string hline =
    "*-" + std::string(max_colsize_name, '-') + "-*-" + std::string(max_colsize_val, '-') + "-*";
  const std::string name_header = "| Post-Processor Name" +
                                  std::string(max_colsize_name - 19, ' ') + " | Value" +
                                  std::string(max_colsize_val - 5, ' ') + " |";

  outstr << hline << "\n";
  outstr << name_header << "\n";
  outstr << hline << "\n";
  for (size_t p = 0; p < num_pps; ++p)
  {
    const auto& [name, val] = scalar_ppnames_and_vals[p];
    outstr << "| " << name << "(latest)" << std::string(max_colsize_name - name.size() - 8, ' ');
    outstr << " | " << std::string(max_colsize_val - val.size(), ' ') << val << " |\n";
  } // for p
  outstr << hline << "\n";

  return outstr.str();
}

void
PostProcessorPrinter::PrintPPsTimeHistory(const std::string& pps_typename,
                                          const std::vector<const PostProcessor*>& pp_list,
                                          const Event& event,
                                          bool per_column_sizes) const
{
  if (pp_list.empty()) return;
  // Establish unique time history sizes
  std::set<size_t> unq_time_histsizes;

  for (const auto& pp : pp_list)
  {
    const size_t time_histsize = pp->GetTimeHistory().size();
    unq_time_histsizes.insert(time_histsize);
  }

  // Subscribe pps to unique time hist sizes
  std::map<size_t, std::vector<const PostProcessor*>> pp_timehist_size_subs;
  for (size_t time_histsize : unq_time_histsizes)
  {
    auto& subs = pp_timehist_size_subs[time_histsize];
    for (const auto& pp : pp_list)
      if (pp->GetTimeHistory().size() == time_histsize) subs.push_back(pp);
  }

  // For each timeline. Build the table
  for (const auto& [timehistsize, pp_sub_list] : pp_timehist_size_subs)
  {
    if (pp_sub_list.empty()) continue;

    //+2 top header + bottom header
    const size_t num_rows = std::min(size_t(time_history_limit_ + 2), timehistsize + 2);
    const size_t num_cols = pp_sub_list.size() + 1; //+1 time column

    auto value_matrix = BuildPPHistoryMatrix(timehistsize, time_history_limit_, pp_sub_list);

    // Find largest column
    size_t max_column_width = 0;
    for (const auto& row : value_matrix)
      for (const auto& entry : row)
        max_column_width = std::max(max_column_width, entry.size());
    max_column_width = std::max(max_column_width, size_t(15)); // minimum size 15

    std::vector<size_t> col_sizes;
    if (not per_column_sizes) col_sizes.assign(num_cols, max_column_width);
    else
    {
      col_sizes.assign(num_cols, 0);
      for (size_t j = 0; j < num_cols; ++j)
      {
        col_sizes[j] = value_matrix[0][j].size();
        for (size_t t = 1; t < num_rows; ++t)
          col_sizes[j] = std::max(col_sizes[j], value_matrix[t][j].size());
      }
    }

    /**Lambda to left pad an entry.*/
    auto LeftPad = [](std::string& entry, size_t width)
    {
      const size_t pad_size = width - entry.size();
      entry = std::string(pad_size, ' ').append(entry);
    };
    /**Lambda to right pad an entry.*/
    auto RightPad = [](std::string& entry, size_t width)
    {
      const size_t pad_size = width - entry.size();
      entry.append(std::string(pad_size, ' '));
    };

    for (size_t j = 0; j < num_cols; ++j)
      RightPad(value_matrix[0][j], col_sizes[j]);
    for (size_t i = 1; i < num_rows; ++i)
      for (size_t j = 0; j < num_cols; ++j)
        LeftPad(value_matrix[i][j], col_sizes[j]);

    // Build sub matrix structure
    std::vector<size_t> sub_mat_sizes;
    size_t total_width = 0;
    size_t col_counter = 0;
    for (size_t c = 0; c < num_cols; ++c)
    {
      //[0]  |           0.000000 |           1.000000 |       0.000000e+00 |
      // 5 chars for log, 2 for '| ', 1 for ' ', 2 for ' |' at end
      const size_t projected_total_width = total_width + col_sizes[c] + 5 + 2 + 1 + 2;
      if (projected_total_width > table_column_limit_)
      {
        total_width = col_sizes[c] + 5 + 2 + 1; // the time column
        sub_mat_sizes.push_back(col_counter);
      }
      col_counter += 1;
      total_width += col_sizes[c] + 2 + 1;
    }
    sub_mat_sizes.push_back(col_counter);

    // Now actually build the sub-matrices
    typedef std::vector<std::string> VecStr;
    typedef std::vector<VecStr> MatStr;
    std::vector<MatStr> sub_matrices;
    size_t last_k = 1;
    for (const size_t k : sub_mat_sizes)
    {
      MatStr sub_matrix(num_rows, VecStr(k - last_k + 1, ""));
      // Copy time col
      for (size_t i = 0; i < num_rows; ++i)
        sub_matrix[i][0] = value_matrix[i][0];

      // Copy other
      for (size_t i = 0; i < num_rows; ++i)
      {
        size_t j_star = 1;
        for (size_t j = last_k; j < k; ++j, ++j_star)
          sub_matrix[i][j_star] = value_matrix[i][j];
      }
      last_k = k;
      sub_matrices.push_back(std::move(sub_matrix));
    }

    std::stringstream outstr;
    for (const auto& sub_matrix : sub_matrices)
      outstr << PrintPPsSubTimeHistory(sub_matrix);

    log.Log() << "\n"
              << pps_typename << " post-processors history at event \"" << event.Name() << "\"\n"
              << outstr.str();
  } // for each thing in pp_timehist_size_subs
}

std::string
PostProcessorPrinter::PrintPPsSubTimeHistory(
  const std::vector<std::vector<std::string>>& sub_history)
{
  const size_t num_rows = sub_history.size();
  const size_t num_cols = sub_history.front().size();

  std::stringstream output;
  std::stringstream hline;
  for (size_t k = 0; k < num_cols; ++k)
  {
    const size_t col_str_size = sub_history.front()[k].size();
    hline << "*" << std::string(col_str_size + 2, '-');
  }
  hline << "*\n";

  output << hline.str();
  for (size_t i = 0; i < num_rows; ++i)
  {
    if (i == 1) output << hline.str();
    std::stringstream line;
    for (size_t k = 0; k < num_cols; ++k)
      line << "| " << sub_history[i][k] << " ";
    line << "|\n";
    output << line.str();
    if (i == (num_rows - 1)) output << hline.str();
  }

  return output.str();
}

void
PostProcessorPrinter::PrintCSVFile(const Event& event) const
{
  const auto scalar_pps = GetScalarPostProcessorsList(event);
  const auto vector_pps = GetVectorPostProcessorsList(event);
  const auto arbitr_pps = GetArbitraryPostProcessorsList(event);

  std::ofstream csvfile;
  csvfile.open(csv_filename_, std::ios::out);

  PrintScalarPPsToCSV(csvfile, scalar_pps);
  PrintVectorPPsToCSV(csvfile, vector_pps);
  PrintArbitraryPPsToCSV(csvfile, arbitr_pps);

  csvfile.close();
}

void
PostProcessorPrinter::PrintScalarPPsToCSV(std::ofstream& csvfile,
                                          const std::vector<const PostProcessor*>& pp_list)
{
  csvfile << "Scalar Post-Processors\n";

  // Establish unique time history sizes
  std::set<size_t> unq_time_histsizes;

  for (const auto& pp : pp_list)
  {
    const size_t time_histsize = pp->GetTimeHistory().size();
    unq_time_histsizes.insert(time_histsize);
  }

  // Subscribe pps to unique time hist sizes
  std::map<size_t, std::vector<const PostProcessor*>> pp_timehist_size_subs;
  for (size_t time_histsize : unq_time_histsizes)
  {
    auto& subs = pp_timehist_size_subs[time_histsize];
    for (const auto& pp : pp_list)
      if (pp->GetTimeHistory().size() == time_histsize) subs.push_back(pp);
  }

  // For each timeline. Build the table
  for (const auto& [timehistsize, pp_sub_list] : pp_timehist_size_subs)
  {
    const auto value_matrix = BuildPPHistoryMatrix(timehistsize, timehistsize, pp_sub_list);
    for (const auto& row : value_matrix)
    {
      for (const auto& entry : row)
      {
        csvfile << entry;
        if (&entry != &row.back()) csvfile << ",";
      }
      csvfile << "\n";
    }
  } // for each thing in pp_timehist_size_subs
}

void
PostProcessorPrinter::PrintVectorPPsToCSV(std::ofstream& csvfile,
                                          const std::vector<const PostProcessor*>& pp_list)
{
  csvfile << "Vector Post-Processors\n";

  for (const auto& pp : pp_list)
  {
    csvfile << pp->Name() << "\n";
    const size_t timehistsize = pp->GetTimeHistory().size();
    const auto value_matrix = BuildPPHistoryMatrix(timehistsize, timehistsize, {pp});
    for (const auto& row : value_matrix)
    {
      for (const auto& entry : row)
      {
        auto entry_star = entry;
        std::replace(entry_star.begin(), entry_star.end(), ' ', ',');
        csvfile << entry_star;
        if (&entry != &row.back()) csvfile << ",";
      }
      csvfile << "\n";
    }
  }
}

void
PostProcessorPrinter::PrintArbitraryPPsToCSV(std::ofstream& csvfile,
                                             const std::vector<const PostProcessor*>& pp_list)
{
  csvfile << "Arbitrary Post-Processors\n";

  for (const auto& pp : pp_list)
  {
    csvfile << pp->Name() << "\n";
    const size_t timehistsize = pp->GetTimeHistory().size();
    const auto value_matrix = BuildPPHistoryMatrix(timehistsize, timehistsize, {pp});
    for (const auto& row : value_matrix)
    {
      for (const auto& entry : row)
      {
        auto entry_star = entry;

        csvfile << entry_star;
        if (&entry != &row.back()) csvfile << ",";
      }
      csvfile << "\n";
    }
  }
}

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

} // namespace opensn
