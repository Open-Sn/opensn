#include "PostProcessorPrinter.h"

#include "post_processors/PostProcessor.h"
#include "event_system/Event.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi
{
// ##################################################################
void PostProcessorPrinter::PrintPPsLatestValuesOnly(
  const std::string& pps_typename,
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
    Chi::log.Log() << "\n"
                   << pps_typename
                   << " post-processors latest values at event \""
                   << event.Name() << "\"\n"
                   << outstr.str() << "\n";
  }
}

// ##################################################################
std::string PostProcessorPrinter::PrintPPsHorizontal(
  const std::vector<std::pair<std::string, std::string>>&
    scalar_ppnames_and_vals,
  int)
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

// ##################################################################
std::string PostProcessorPrinter::PrintPPsVertical(
  const std::vector<std::pair<std::string, std::string>>&
    scalar_ppnames_and_vals,
  int event_code)
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

  const std::string hline = "*-" + std::string(max_colsize_name, '-') + "-*-" +
                            std::string(max_colsize_val, '-') + "-*";
  const std::string name_header =
    "| Post-Processor Name" + std::string(max_colsize_name - 19, ' ') +
    " | Value" + std::string(max_colsize_val - 5, ' ') + " |";

  outstr << hline << "\n";
  outstr << name_header << "\n";
  outstr << hline << "\n";
  for (size_t p = 0; p < num_pps; ++p)
  {
    const auto& [name, val] = scalar_ppnames_and_vals[p];
    outstr << "| " << name << "(latest)"
           << std::string(max_colsize_name - name.size() - 8, ' ');
    outstr << " | " << std::string(max_colsize_val - val.size(), ' ') << val
           << " |\n";
  } // for p
  outstr << hline << "\n";

  return outstr.str();
}
} // namespace chi
