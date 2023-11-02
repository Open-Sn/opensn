#include "PostProcessorPrinter.h"

#include "post_processors/PostProcessor.h"
#include "chi_utils.h"

#include <set>
#include <algorithm>

namespace chi
{

void PostProcessorPrinter::PrintCSVFile(const Event& event) const
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

void PostProcessorPrinter::PrintScalarPPsToCSV(
  std::ofstream& csvfile, const std::vector<const PostProcessor*>& pp_list)
{
  csvfile << "Scalar Post-Processors\n";

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
    const auto value_matrix =
      BuildPPHistoryMatrix(timehistsize, timehistsize, pp_sub_list);
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

void PostProcessorPrinter::PrintVectorPPsToCSV(
  std::ofstream& csvfile, const std::vector<const PostProcessor*>& pp_list)
{
  csvfile << "Vector Post-Processors\n";

  for (const auto& pp : pp_list)
  {
    csvfile << pp->Name() << "\n";
    const size_t timehistsize = pp->GetTimeHistory().size();
    const auto value_matrix =
      BuildPPHistoryMatrix(timehistsize, timehistsize, {pp});
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

void PostProcessorPrinter::PrintArbitraryPPsToCSV(
  std::ofstream& csvfile, const std::vector<const PostProcessor*>& pp_list)
{
  csvfile << "Arbitrary Post-Processors\n";

  for (const auto& pp : pp_list)
  {
    csvfile << pp->Name() << "\n";
    const size_t timehistsize = pp->GetTimeHistory().size();
    const auto value_matrix =
      BuildPPHistoryMatrix(timehistsize, timehistsize, {pp});
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

} // namespace chi