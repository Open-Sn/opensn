#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"

#include "framework/logging/stringstream_color.h"

#include <sstream>

namespace opensn
{

Logger&
Logger::GetInstance() noexcept
{
  static Logger instance;
  return instance;
}

Logger::Logger() noexcept
{
  verbosity_ = 0;
  std::string memory_usage_event("Maximum Memory Usage");
  repeating_events.emplace_back(memory_usage_event);

  RepeatingEvent& ref_rep_event = repeating_events.back();

  ref_rep_event.Events().emplace_back(
    program_timer.GetTime(), EventType::EVENT_CREATED, std::make_shared<EventInfo>());
}

LogStream
Logger::Log(LOG_LVL level)
{
  switch (level)
  {
    case LOG_0VERBOSE_0:
    case LOG_0:
    {
      if (opensn::mpi_comm.rank() == 0)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_0WARNING:
    {
      if (opensn::mpi_comm.rank() == 0)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_YELLOW) + "**WARNING** ";
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_0ERROR:
    {
      if (opensn::mpi_comm.rank() == 0)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_RED) + "**!**ERROR**!** ";
        return {&std::cerr, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }

    case LOG_0VERBOSE_1:
    {
      if ((opensn::mpi_comm.rank() == 0) and (verbosity_ >= 1))
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_CYAN);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case Logger::LOG_LVL::LOG_0VERBOSE_2:
    {
      if ((opensn::mpi_comm.rank() == 0) and (verbosity_ >= 2))
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_MAGENTA);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_ALLVERBOSE_0:
    case LOG_ALL:
    {
      std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
      return {&std::cout, header};
    }
    case LOG_ALLWARNING:
    {
      std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
      header += StringStreamColor(FG_YELLOW) + "**WARNING** ";
      return {&std::cout, header};
    }
    case LOG_ALLERROR:
    {
      std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
      header += StringStreamColor(FG_RED) + "**!**ERROR**!** ";
      return {&std::cerr, header};
    }

    case LOG_ALLVERBOSE_1:
    {
      if (verbosity_ >= 1)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_CYAN);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    case LOG_ALLVERBOSE_2:
    {
      if (verbosity_ >= 2)
      {
        std::string header = "[" + std::to_string(opensn::mpi_comm.rank()) + "]  ";
        header += StringStreamColor(FG_MAGENTA);
        return {&std::cout, header};
      }
      else
      {
        std::string header = " ";
        return {&dummy_stream_, header, true};
      }
    }
    default:
      std::string header = " ";
      return {&dummy_stream_, header};
  }
}

void
Logger::SetVerbosity(int int_level)
{
  verbosity_ = std::min(int_level, 2);
}

int
Logger::GetVerbosity() const
{
  return verbosity_;
}

size_t
Logger::GetRepeatingEventTag(std::string event_name)
{
  repeating_events.emplace_back(event_name);

  RepeatingEvent& ref_rep_event = repeating_events.back();

  ref_rep_event.Events().emplace_back(
    program_timer.GetTime(), EventType::EVENT_CREATED, std::make_shared<EventInfo>());

  return repeating_events.size() - 1;
}

size_t
Logger::GetExistingRepeatingEventTag(std::string event_name)
{
  const size_t num_rep_events = repeating_events.size();
  for (size_t k = num_rep_events - 1; k != 0; --k)
    if (repeating_events[k].Name() == event_name)
      return k;

  OpenSnLogicalError("Tag could not be found for repeating event name \"" + event_name + "\"");
}

void
Logger::LogEvent(size_t ev_tag, EventType ev_type, const std::shared_ptr<EventInfo>& ev_info)
{
  if (ev_tag >= repeating_events.size())
    return;

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  ref_rep_event.Events().emplace_back(program_timer.GetTime(), ev_type, ev_info);
}

void
Logger::LogEvent(size_t ev_tag, EventType ev_type)
{
  if (ev_tag >= repeating_events.size())
    return;

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  ref_rep_event.Events().emplace_back(program_timer.GetTime(), ev_type, nullptr);
}

std::string
Logger::PrintEventHistory(size_t ev_tag)
{
  std::stringstream outstr;
  if (ev_tag >= repeating_events.size())
    return outstr.str();

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  for (auto& event : ref_rep_event.Events())
  {
    outstr << "[" << opensn::mpi_comm.rank() << "] ";

    char buf[100];
    snprintf(buf, 100, "%16.9f", event.ev_time / 1000.0);
    outstr << buf << " ";

    switch (event.ev_type)
    {
      case EventType::EVENT_CREATED:
        outstr << "EVENT_CREATED ";
        break;
      case EventType::SINGLE_OCCURRENCE:
        outstr << "SINGLE_OCCURRENCE ";
        break;
      case EventType::EVENT_BEGIN:
        outstr << "EVENT_BEGIN ";
        break;
      case EventType::EVENT_END:
        outstr << "EVENT_END ";
        break;
    }

    if (event.ev_info != nullptr)
      outstr << event.ev_info->GetString();
    outstr << std::endl;
  }

  return outstr.str();
}

double
Logger::ProcessEvent(size_t ev_tag, Logger::EventOperation ev_operation)
{
  if (ev_tag >= repeating_events.size())
    return 0.0;

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  double ret_val = 0.0;
  switch (ev_operation)
  {
    case EventOperation::NUMBER_OF_OCCURRENCES:
    {
      for (auto& event : ref_rep_event.Events())
      {
        if ((event.ev_type == EventType::EVENT_CREATED) or
            (event.ev_type == EventType::SINGLE_OCCURRENCE) or
            (event.ev_type == EventType::EVENT_BEGIN))
          ret_val += 1.0;
      } // for events
      break;
    }
    case EventOperation::TOTAL_DURATION:
    {
      double start_time = 0.0;
      for (auto& event : ref_rep_event.Events())
      {
        if (event.ev_type == EventType::EVENT_BEGIN)
          start_time = event.ev_time;
        if (event.ev_type == EventType::EVENT_END)
          ret_val += event.ev_time - start_time;
      } // for events
      ret_val *= 1000.0;
      break;
    }
    case EventOperation::AVERAGE_DURATION:
    {
      double start_time = 0.0;
      int counter = 0;
      for (auto& event : ref_rep_event.Events())
      {
        if (event.ev_type == EventType::EVENT_BEGIN)
          start_time = event.ev_time;

        if (event.ev_type == EventType::EVENT_END)
        {
          ret_val += event.ev_time - start_time;
          counter++;
        }
      } // for events
      ret_val /= (1000.0 * counter);
      break;
    }
    case EventOperation::MAX_VALUE:
    {
      ret_val = 0.0;
      for (auto& event : ref_rep_event.Events())
      {
        if ((event.ev_type == EventType::SINGLE_OCCURRENCE) or
            (event.ev_type == EventType::EVENT_BEGIN) or (event.ev_type == EventType::EVENT_END))
        {
          if (event.ev_info != nullptr)
            ret_val = std::max(event.ev_info->arb_value, ret_val);
        }
      } // for events
      break;
    }
    case EventOperation::AVERAGE_VALUE:
    {
      ret_val = 0.0;
      int count = 0;
      for (auto& event : ref_rep_event.Events())
      {
        if ((event.ev_type == EventType::SINGLE_OCCURRENCE) or
            (event.ev_type == EventType::EVENT_BEGIN) or (event.ev_type == EventType::EVENT_END))
        {
          if (event.ev_info != nullptr)
          {
            ret_val += event.ev_info->arb_value;
            ++count;
          }
        }
      } // for events
      if (count == 0)
        count = 1;
      ret_val /= count;
      break;
    }
  } // switch

  return ret_val;
}

} // namespace opensn
