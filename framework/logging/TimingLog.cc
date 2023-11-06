#include "framework/logging/TimingLog.h"

#include "framework/utils/chi_timer.h"

#include "framework/chi_runtime.h"

#include "framework/logging/chi_log_exceptions.h"

#include <sstream>

namespace chi
{

TimingBlock&
TimingLog::CreateTimingBlock(const std::string& name, const std::string& parent_name /*=""*/)
{
  ChiInvalidArgumentIf(timing_blocks_.count(name) != 0,
                       "TimingBlock with name \"" + name +
                         "\" already exists"
                         " in the logger");
  auto block_ptr = std::make_unique<TimingBlock>(name);
  auto saved_pointer = &(*block_ptr);
  timing_blocks_[name] = std::move(block_ptr);

  if (parent_name.empty())
  {
    if (name != "ChiTech")
    {
      auto iter = timing_blocks_.find("ChiTech");
      ChiLogicalErrorIf(iter == timing_blocks_.end(),
                        "Bad error, could not fine the \"ChiTech\" timing block");

      iter->second->AddChild(*saved_pointer);
    }

    return *saved_pointer;
  }

  auto iter = timing_blocks_.find(parent_name);

  ChiInvalidArgumentIf(iter == timing_blocks_.end(),
                       "Parent name \"" + parent_name + "\" does not exist.");

  auto& parent = iter->second;
  parent->AddChild(*saved_pointer);

  return *saved_pointer;
}

TimingBlock&
TimingLog::CreateOrGetTimingBlock(const std::string& name, const std::string& parent_name /*=""*/)
{
  auto iter = timing_blocks_.find(name);

  if (iter == timing_blocks_.end()) return CreateTimingBlock(name, parent_name);

  return *iter->second;
}

TimingBlock&
TimingLog::GetTimingBlock(const std::string& name)
{
  auto iter = timing_blocks_.find(name);

  ChiInvalidArgumentIf(iter == timing_blocks_.end(),
                       "Timing block with name \"" + name + "\" does not exist.");

  return *iter->second;
}

// ##################################################################

TimingBlock::TimingBlock(const std::string& name) : name_(name)
{
}

void
TimingBlock::TimeSectionBegin()
{
  reference_time_ = Chi::program_timer.GetTime();
}

void
TimingBlock::TimeSectionEnd()
{
  const double delta_t = Chi::program_timer.GetTime() - reference_time_;
  total_time_ += delta_t;
  num_occurences_ += 1;
  last_delta_time_ = delta_t;
}

size_t
TimingBlock::NumberOfOccurences() const
{
  return num_occurences_;
}

double
TimingBlock::TotalTime() const
{
  return total_time_;
}

double
TimingBlock::AverageTime() const
{
  if (num_occurences_ == 0) return 0.0;

  return total_time_ / static_cast<double>(num_occurences_);
}

double
TimingBlock::LastDelta() const
{
  return last_delta_time_;
}

void
TimingBlock::AddChild(const TimingBlock& child_block)
{
  children_.push_back(&child_block);
}

std::string
TimingBlock::MakeGraphString()
{
  if (name_ == "ChiTech") TimeSectionEnd();

  std::vector<std::vector<std::string>> string_matrix;

  AppendGraphEntry(string_matrix, nullptr, "");

  const size_t J = 5;
  const size_t I = string_matrix.size();

  std::vector<size_t> max_col_widths(J, 0);
  for (size_t i = 0; i < I; ++i)
    for (size_t j = 0; j < J; ++j)
      max_col_widths[j] = std::max(max_col_widths[j], string_matrix[i][j].size());

  std::vector<std::string> headers = {
    "Section Name", "#calls", "Total time[s]", "Average time[s]", "% of parent"};

  for (size_t j = 0; j < J; ++j)
    max_col_widths[j] = std::max(max_col_widths[j], headers[j].size());

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

  std::stringstream outstr;

  auto HDIV = [&outstr, &max_col_widths]()
  {
    outstr << "*-";
    for (size_t j = 0; j < J; ++j)
    {
      outstr << std::string(max_col_widths[j] + 1, '-');
      if (j < (J - 1)) outstr << "*-";
    }
    outstr << "*\n";
  };

  HDIV();

  outstr << "| ";
  for (size_t j = 0; j < J; ++j)
  {
    if (j == 0) RightPad(headers[j], max_col_widths[j]);
    else
      LeftPad(headers[j], max_col_widths[j] + 1);
    outstr << headers[j] << " |";
  }
  outstr << "\n";

  HDIV();

  for (size_t i = 0; i < I; ++i)
  {
    outstr << "| ";
    for (size_t j = 0; j < J; ++j)
    {
      if (j == 0) RightPad(string_matrix[i][j], max_col_widths[j]);
      else
        LeftPad(string_matrix[i][j], max_col_widths[j] + 1);
      outstr << string_matrix[i][j] << " |";
    }
    outstr << "\n";
  }

  HDIV();

  return outstr.str();
}

//  NOLINTBEGIN(misc-no-recursion)
void
TimingBlock::AppendGraphEntry(std::vector<std::vector<std::string>>& string_matrix,
                              const TimingBlock* parent,
                              const std::string& indent) const
{
  std::vector<std::string> entry;

  const double parent_total_time = parent ? parent->TotalTime() : 0.0;
  const double relative_time =
    parent_total_time > 1.0e-12 ? total_time_ * 100.0 / parent_total_time : 0.0;

  entry.push_back(indent + name_);
  entry.push_back(std::to_string(num_occurences_));

  {
    char buffer[20];
    ChiLogicalErrorIf(snprintf(buffer, 20, "%.5g", total_time_ / 1000.0) < 0,
                      "Failed to convert total_time = " + std::to_string(total_time_));
    entry.push_back(buffer);
  }

  {
    char buffer[20];
    ChiLogicalErrorIf(snprintf(buffer, 20, "%.5g", AverageTime() / 1000.0) < 0,
                      "Failed to convert AverageTime = " + std::to_string(AverageTime()));
    entry.push_back(buffer);
  }

  {
    char buffer[10];
    ChiLogicalErrorIf(snprintf(buffer, 10, "%.2f%%", relative_time) < 0,
                      "Failed to convert relative_time = " + std::to_string(relative_time));
    entry.push_back(parent ? buffer : "--");
  }

  string_matrix.push_back(std::move(entry));

  for (auto& child_ptr : children_)
    child_ptr->AppendGraphEntry(string_matrix, this, indent + "  ");
}
//  NOLINTEND(misc-no-recursion)

} // namespace chi
