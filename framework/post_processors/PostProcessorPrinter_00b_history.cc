#include "PostProcessorPrinter.h"

#include "event_system/Event.h"

#include "post_processors/PostProcessor.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi
{

// ##################################################################
void PostProcessorPrinter::PrintPPsTimeHistory(
  const std::string& pps_typename,
  const std::vector<const PostProcessor*>& pp_list,
  const Event& event,
  bool per_column_sizes /*=false*/) const
{
  if (pp_list.empty()) return;
  //======================================== Establish unique time history sizes
  std::set<size_t> unq_time_histsizes;

  for (const auto& pp : pp_list)
  {
    const size_t time_histsize = pp->GetTimeHistory().size();
    unq_time_histsizes.insert(time_histsize);
  }

  //======================================== Subscribe pps to unique time
  //                                         hist sizes
  std::map<size_t, std::vector<const PostProcessor*>> pp_timehist_size_subs;
  for (size_t time_histsize : unq_time_histsizes)
  {
    auto& subs = pp_timehist_size_subs[time_histsize];
    for (const auto& pp : pp_list)
      if (pp->GetTimeHistory().size() == time_histsize) subs.push_back(pp);
  }

  //======================================== For each timeline. Build the table
  for (const auto& [timehistsize, pp_sub_list] : pp_timehist_size_subs)
  {
    if (pp_sub_list.empty()) continue;

    //+2 top header + bottom header
    const size_t num_rows =
      std::min(size_t(time_history_limit_ + 2), timehistsize + 2);
    const size_t num_cols = pp_sub_list.size() + 1; //+1 time column

    auto value_matrix =
      BuildPPHistoryMatrix(timehistsize, time_history_limit_, pp_sub_list);

    // Find largest column
    size_t max_column_width = 0;
    for (const auto& row : value_matrix)
      for (const auto& entry : row)
        max_column_width = std::max(max_column_width, entry.size());
    max_column_width =
      std::max(max_column_width, size_t(15)); // minimum size 15

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
      const size_t projected_total_width =
        total_width + col_sizes[c] + 5 + 2 + 1 + 2;
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

    Chi::log.Log() << "\n"
                   << pps_typename << " post-processors history at event \""
                   << event.Name() << "\"\n"
                   << outstr.str();
  } // for each thing in pp_timehist_size_subs
}

std::string PostProcessorPrinter::PrintPPsSubTimeHistory(
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

} // namespace chi