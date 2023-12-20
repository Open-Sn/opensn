#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include <memory>

namespace opensn
{

class TimingBlock;

/**Utility class for defining time logs.*/
class TimingLog
{
public:
  /**Creates a time block and returns a reference to it. If the name
   * is already taken or the parent name is not found then this method
   * throws `std::invalid_argument`. The parent name may be empty, i.e. "",
   * which will revert to the main timing block "ChiTech"*/
  TimingBlock& CreateTimingBlock(const std::string& name, const std::string& parent_name = "");
  /**If the timing block with the given exists, this method returns that block,
   * otherwise it creates a time block and returns a reference to it. If the
   * parent name is not found then this method throws `std::invalid_argument`.
   * The parent name may be empty, i.e. "", which will revert to the main timing
   * block "ChiTech"*/
  TimingBlock& CreateOrGetTimingBlock(const std::string& name, const std::string& parent_name = "");
  /**Returns a reference to the timing block with the given name. If a timing
   * block with that name does not exist then this method will throw
   * `std::invalid_argument`.*/
  TimingBlock& GetTimingBlock(const std::string& name);

protected:
  std::map<std::string, std::unique_ptr<TimingBlock>> timing_blocks_;
};

/**Hierarchical timing block to efficiently capture timing.*/
class TimingBlock
{
public:
  /**Constructs the given block with the given name.*/
  explicit TimingBlock(const std::string& name);
  /**Deleted copy constructor.*/
  TimingBlock(const TimingBlock& other) = delete;
  /**Deleted move constructor.*/
  TimingBlock(TimingBlock&& other) = delete;

  /**Begins a timing section by resetting the reference time.*/
  void TimeSectionBegin();
  /**Ends a timing section, contributes to the total time and adds an
   * occurrence.*/
  void TimeSectionEnd();

  /**Returns the number of timing section occurrences.*/
  size_t NumberOfOccurences() const;
  /**Returns the computed total time over all occurrences.*/
  double TotalTime() const;
  /**Returns the average time per occurence.*/
  double AverageTime() const;
  /**Returns the last computed time differential.*/
  double LastDelta() const;

  /**Makes a string table of the timing block and all its children.*/
  std::string MakeGraphString();

protected:
  friend class TimingLog;
  /**Adds the supplied timing black as a child.*/
  void AddChild(const TimingBlock& child_block);
  /**Used when building the graph, this is a recursive function that adds
   * each entry to a matrix.*/
  void AppendGraphEntry(std::vector<std::vector<std::string>>& string_matrix,
                        const TimingBlock* parent,
                        const std::string& indent) const;

  const std::string name_;
  size_t num_occurences_ = 0;
  double total_time_ = 0.0;
  double reference_time_ = 0.0;
  double last_delta_time_ = 0.0;
  std::vector<const TimingBlock*> children_;
};

} // namespace opensn
