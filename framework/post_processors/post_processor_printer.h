#pragma once

#include "framework/event_system/event_subscriber.h"
#include <memory>
#include <vector>
#include <string>
#include <fstream>

namespace chi
{
class PostProcessor;
class PostProcessorPrinter;
class ParameterBlock;

/**
 * A helper object to allow the printer to subscribe to events.
 */
class PPPrinterSubscribeHelper : public EventSubscriber
{
public:
  explicit PPPrinterSubscribeHelper(PostProcessorPrinter& printer_ref);

  void ReceiveEventUpdate(const Event& event) override;

private:
  PostProcessorPrinter& printer_ref_;
};

enum class ScalarPPTableFormat : int
{
  VERTICAL = 0,
  HORIZONTAL = 1,
};

/**
 * A singleton responsible for printing post-processors.
 */
class PostProcessorPrinter
{
public:
  static PostProcessorPrinter& GetInstance();

  PostProcessorPrinter(const PostProcessorPrinter&) = delete;
  PostProcessorPrinter operator=(const PostProcessorPrinter&) = delete;

  void ReceiveEventUpdate(const Event& event);

  static char SubscribeToSystemWideEventPublisher();

  void SetScalarPPTableFormat(ScalarPPTableFormat format);
  void SetEventsOnWhichPrintPPs(const std::vector<std::string>& events);

  void SetPrintScalarTimeHistory(bool value);
  void SetPrintVectorTimeHistory(bool value);

  void SetScalarPerColumnSize(bool value);
  void SetVectorPerColumnSize(bool value);

  void SetTableColumnLimit(size_t limit);
  void SetTimeHistoryLimit(size_t limit);

  void SetCSVFilename(const std::string& csv_filename);

  /**
   * A manual means to print a post processor.
   */
  std::string GetPrintedPostProcessors(const std::vector<const PostProcessor*>& pp_list) const;

private:
  PostProcessorPrinter();

  void PrintPostProcessors(const Event& event) const;

  void PrintPPsLatestValuesOnly(const std::string& pps_typename,
                                const std::vector<const PostProcessor*>& pp_list,
                                const Event& event) const;

  static std::string PrintPPsHorizontal(
    const std::vector<std::pair<std::string, std::string>>& scalar_ppnames_and_vals, int);
  static std::string
  PrintPPsVertical(const std::vector<std::pair<std::string, std::string>>& scalar_ppnames_and_vals,
                   int event_code);

  void PrintPPsTimeHistory(const std::string& pps_typename,
                           const std::vector<const PostProcessor*>& pp_list,
                           const Event& event,
                           bool per_column_sizes = false) const;

  static std::string
  PrintPPsSubTimeHistory(const std::vector<std::vector<std::string>>& sub_history);

  void PrintCSVFile(const Event& event) const;
  static void PrintScalarPPsToCSV(std::ofstream& csvfile,
                                  const std::vector<const PostProcessor*>& pp_list);
  static void PrintVectorPPsToCSV(std::ofstream& csvfile,
                                  const std::vector<const PostProcessor*>& pp_list);
  static void PrintArbitraryPPsToCSV(std::ofstream& csvfile,
                                     const std::vector<const PostProcessor*>& pp_list);

  static std::vector<const PostProcessor*> GetScalarPostProcessorsList(const Event& event);

  static std::vector<const PostProcessor*> GetVectorPostProcessorsList(const Event& event);

  static std::vector<const PostProcessor*> GetArbitraryPostProcessorsList(const Event& event);

  static std::vector<std::vector<std::string>>
  BuildPPHistoryMatrix(size_t timehistsize,
                       size_t time_history_limit,
                       const std::vector<const PostProcessor*>& pp_sub_list);

  static std::shared_ptr<PPPrinterSubscribeHelper> helper_ptr_;
  std::vector<std::string> events_on_which_to_print_postprocs_;

  ScalarPPTableFormat scalar_pp_table_format_ = ScalarPPTableFormat::VERTICAL;
  bool print_scalar_time_history_ = true;
  bool print_vector_time_history_ = true;
  bool per_column_size_scalars_ = true;
  bool per_column_size_vectors_ = true;
  size_t table_column_limit_ = 120;
  size_t time_history_limit_ = 15;

  std::string csv_filename_;
};

} // namespace chi
