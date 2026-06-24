// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <iostream>
#include <sstream>
#include <string_view>
#include <utility>
#include "stringstream_color.h"

namespace opensn
{

inline int
NoWrapStreamIndex()
{
  static const int index = std::ios_base::xalloc();
  return index;
}

inline std::ostream&
no_wrap(std::ostream& stream)
{
  stream.iword(NoWrapStreamIndex()) = 1;
  return stream;
}

/// Log stream for adding header information to a string stream.
class LogStream : public std::stringstream
{
public:
  struct Header
  {
    std::string prefix;
    std::string label;
    std::string color;
  };

  LogStream(std::ostream* output_stream,
            Header header,
            bool dummy_flag = false,
            bool use_color = false)
    : log_stream_(output_stream),
      header_(std::move(header)),
      dummy_(dummy_flag),
      use_color_(use_color)
  {
  }

  LogStream(std::ostream* output_stream, std::string header, bool dummy_flag = false)
    : LogStream(output_stream, Header{std::move(header), "", ""}, dummy_flag, false)
  {
  }

  LogStream(const LogStream&) = delete;
  LogStream& operator=(const LogStream&) = delete;

  ~LogStream() noexcept override
  {
    if (dummy_)
      return;

    try
    {
      std::string content = this->str();
      if (content.empty())
        return;

      std::istringstream iss(content);
      std::string line;
      std::string oline;
      std::string reset_str = use_color_ ? StringStreamColor(StringStreamColorCode::RESET) : "";
      const bool wrap_lines = this->iword(NoWrapStreamIndex()) == 0;
      while (std::getline(iss, line))
      {
        if (wrap_lines)
          AppendWrappedLine(oline, line, reset_str);
        else
          AppendUnwrappedLine(oline, line, reset_str);
      }

      if (!oline.empty())
        *log_stream_ << oline << std::flush;
    }
    catch (...) // NOLINT(bugprone-empty-catch)
    {
      // No exceptions escape the destructor...
    }
  }

private:
  static constexpr std::size_t LINE_WRAP_WIDTH = 150;
  static constexpr std::size_t CONTINUATION_INDENT = 4;

  static std::size_t FindWrapPosition(std::string_view text, std::size_t max_width)
  {
    if (text.size() <= max_width)
      return std::string_view::npos;

    const auto search_limit = std::min(max_width, text.size() - 1);
    for (std::size_t i = search_limit; i > 0; --i)
      if (text[i] == ' ')
        return i;

    for (std::size_t i = max_width + 1; i < text.size(); ++i)
      if (text[i] == ' ')
        return i;

    return std::string_view::npos;
  }

  static std::string_view TrimLeadingSpaces(std::string_view text)
  {
    while (not text.empty() and text.front() == ' ')
      text.remove_prefix(1);

    return text;
  }

  void AppendUnwrappedLine(std::string& output,
                           std::string_view line,
                           const std::string& reset_str) const
  {
    output += header_.prefix;
    output += header_.color;
    output += header_.label;
    output += line;
    output += reset_str;
    output += "\n";
  }

  void
  AppendWrappedLine(std::string& output, std::string_view line, const std::string& reset_str) const
  {
    bool first_piece = true;

    bool line_remaining = true;
    while (line_remaining)
    {
      const auto header_width = first_piece ? header_.prefix.size() + header_.label.size()
                                            : header_.prefix.size() + CONTINUATION_INDENT;
      const auto content_width =
        header_width < LINE_WRAP_WIDTH ? LINE_WRAP_WIDTH - header_width : LINE_WRAP_WIDTH;
      const auto wrap_pos = FindWrapPosition(line, content_width);
      const auto piece = wrap_pos == std::string_view::npos ? line : line.substr(0, wrap_pos);

      output += header_.prefix;
      output += header_.color;
      if (first_piece)
        output += header_.label;
      else
        output += std::string(CONTINUATION_INDENT, ' ');

      output += piece;
      output += reset_str;
      output += "\n";

      first_piece = false;
      line_remaining = wrap_pos != std::string_view::npos;

      if (line_remaining)
        line = TrimLeadingSpaces(line.substr(wrap_pos));
    }
  }

  std::ostream* log_stream_;
  Header header_;
  const bool dummy_;
  bool use_color_;
};

struct DummyStream : public std::ostream
{
  struct DummyStreamBuffer : std::streambuf
  {
  protected:
    int overflow(int c) override { return c; };
  } buffer;

  DummyStream() : std::ostream(&buffer) {}
  ~DummyStream() override = default;
};

} // namespace opensn
