#include "opensn/framework/logging/chi_logstream.h"

#include "opensn/framework/logging/stringstream_color.h"

chi::LogStream::~LogStream()
{
  if (dummy_) return;

  std::string line, oline;
  while (std::getline(*this, line))
    oline += log_header_ + line + '\n' + StringStreamColor(RESET);

  if (!oline.empty()) *log_stream_ << oline << std::flush;
}
