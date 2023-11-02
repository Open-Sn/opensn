#ifndef CHI_LOGSTREAM_H
#define CHI_LOGSTREAM_H

#include <iostream>
#include <sstream>

namespace chi
{
//###################################################################
/** Log stream for adding header information to a string stream.*/
class LogStream : public std::stringstream
{
private:
  std::ostream* log_stream_;
  std::string  log_header_;
  const bool dummy_ = false;

public:
  /** Creates a string stream.*/
  LogStream(std::ostream* output_stream,
            std::string header,
            bool dummy_flag=false) :
    log_stream_(output_stream),
    log_header_(std::move(header)),
    dummy_(dummy_flag)
  { }

  /** Flushes the broken-up/headered stream to the output.*/
  virtual ~LogStream();

  LogStream(const LogStream& other)
  {
    log_stream_ = other.log_stream_;
    log_header_ = other.log_header_;
  }
};

struct DummyStream: public std::ostream
{
  struct DummyStreamBuffer : std::streambuf
  {
    virtual int overflow(int c) { return c; };
  } buffer;

  DummyStream(): std::ostream(&buffer) {}

  ~DummyStream() {}
};
}//namespace chi_objects
#endif //CHI_LOGSTREAM_H