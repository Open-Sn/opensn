// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <iostream>
#include <sstream>

namespace opensn
{

/// Log stream for adding header information to a string stream.
class LogStream : public std::stringstream
{
private:
  std::ostream* log_stream_;
  std::string log_header_;
  const bool dummy_ = false;

public:
  /// Creates a string stream.
  LogStream(std::ostream* output_stream, std::string header, bool dummy_flag = false)
    : log_stream_(output_stream), log_header_(std::move(header)), dummy_(dummy_flag)
  {
  }

  /// Flushes stream.
  virtual ~LogStream();

  LogStream(const LogStream& other)
  {
    log_stream_ = other.log_stream_;
    log_header_ = other.log_header_;
  }
};

struct DummyStream : public std::ostream
{
  struct DummyStreamBuffer : std::streambuf
  {
    virtual int overflow(int c) { return c; };
  } buffer;

  DummyStream() : std::ostream(&buffer) {}

  ~DummyStream() {}
};

} // namespace opensn
