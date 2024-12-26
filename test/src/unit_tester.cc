#include "test/src/unit_tester.h"

#include "framework/logging/log.h"
#include "framework/mpi/mpi_utils.h"

using namespace opensn;

namespace unit_tests
{

UnitTester::UnitTester() : passed_(0), failed_(0)
{
}

void
UnitTester::test(const std::string& file,
                 const int line,
                 const std::optional<std::string>& failed,
                 const bool throw_on_fail)
{
  if (failed.has_value())
  {
    const std::string message =
      file + ":" + std::to_string(line) + ": FAILED\n  REASON: " + *failed;
    if (throw_on_fail)
      throw std::logic_error(message);
    ++failed_;
    opensn::log.LogAllVerbose0() << message;
  }
  else
  {
    ++passed_;
    opensn::log.LogAllVerbose2() << file << ":" << line << ": PASSED";
  }
}

void
UnitTester::finalize() const
{
  opensn::log.LogAllVerbose0() << "Ran " << passed_ + failed_ << " unit tests; " << failed_
                               << " failed";

  unsigned int total_passed, total_failed;
  opensn::mpi_comm.all_reduce(passed_, total_passed, mpi::op::sum<unsigned int>());
  opensn::mpi_comm.all_reduce(failed_, total_failed, mpi::op::sum<unsigned int>());

  if (opensn::mpi_comm.rank() == 0)
    opensn::log.Log0() << "Ran " << total_passed + total_failed << " total unit tests; "
                       << total_failed << " failed";

  if (total_failed > 0)
    throw std::logic_error("Unit tests failed");
}

} // namespace unit_tests
