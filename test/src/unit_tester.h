#pragma once

#include "framework/runtime.h"

#include "lua/framework/console/console.h"

#include <string>
#include <optional>
#include <type_traits>

using namespace opensn;

/// Sets up the unit tester (call once within a method first)
#define SETUP_UNIT_TESTER() UnitTester tester
/// Finalizes the unit tester which will throw when tests fail (call at the end)
#define FINALIZE_UNIT_TESTER() tester.finalize()
/// Tests if values are equal
#define EXPECT_EQUAL(expected, actual)                                                             \
  tester.test_equal(__FILE__, __LINE__, expected, actual, false)
/// Tests if a value is true
#define EXPECT_TRUE(condition) tester.test_equal(__FILE__, __LINE__, true, condition, false);

namespace unit_tests
{
/**
 * Helper class that can be utilized in a standard LUA test in a more unit-like manner.
 *
 * Usage example:
 *
 *   SETUP_UNIT_TESTER();
 *   EXPECT_EQUAL(false, true); // will fail; false != true
 *   EXPECT_TRUE(true); // will pass
 *   FINALIZE_UNIT_TESTER(); // will throw; 1 test failed
 *
 * It will output file location and line as errors are hit. As errors are collected,
 * it will allow testing to continue (instead of exiting on a failure) and will
 * throw upon finalize if errors are found (also in parallel). For example:
 *
 * ...
 * [1]  opensn/test/framework/mesh/point_inside_cell/point_inside_cell_test.cc:32: FAILED
 * [1]    REASON: Expected value 1 != 0
 * [0]  Ran 162 unit tests; 18 failed
 *
 * The methods on this class should not be used directly. The macros should be utilized
 * as they allow for the collection of fails with additional context (file and line).
 */
class UnitTester
{
public:
  UnitTester();

  /**
   * Internal method for testing equality. Should be used via the EXPECT_EQUAL() macro.
   *
   * Will output the value to screen if T supports the ostream << overload.
   *
   * @param file Path to the file (obtained via a macro)
   * @param line Line in the file (obtained via a macro)
   * @param expected The expected value
   * @param actual The actual value
   * @param bool throw_on_fail Whether or not to throw during a failure
   */
  template <typename T>
  void test_equal(const std::string& file,
                  const int line,
                  const T& expected,
                  const T& actual,
                  const bool throw_on_fail);

  /// Internal method for finalizing. Errors if any tests have failed.
  /// Should be used via the FINALIZE_UNIT_TESTER() macro.
  void finalize() const;

  /// A helper trait to determine if a type supports streaming. This
  /// enables us to output data types to screen if possible during a failure.
  template <class T>
  class is_streamable
  {
    template <class TT>
    static auto test(int)
      -> decltype(std::declval<std::ostream&>() << std::declval<TT>(), std::true_type());

    template <class>
    static auto test(...) -> std::false_type;

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
  };

  /// Whether or not the given type supports streaming, i.e.:
  /// is_streamable_v<double> = true (you can do "<< value" for a double)
  template <class T>
  inline static constexpr bool is_streamable_v = is_streamable<T>::value;

private:
  /**
   * Internal method that tests use to capture results.
   * @param file Path to the file (obtained via a macro)
   * @param line Line in the file (obtained via a macro)
   * @param failed The reason the test failed (if any)
   * @param throw_on_fail Whether or not to throw during a failure
   */
  void test(const std::string& file,
            const int line,
            const std::optional<std::string>& failed,
            const bool throw_on_fail);

  /// Number of tests that have passed locally
  unsigned int passed_;
  /// Number of tests that have failed locally
  unsigned int failed_;
};

template <typename T>
void
UnitTester::test_equal(const std::string& file,
                       const int line,
                       const T& expected,
                       const T& actual,
                       const bool throw_on_fail)
{
  std::optional<std::string> failed;
  if (expected != actual)
  {
    if constexpr (is_streamable_v<T>) // output value if we can!
    {
      std::stringstream reason;
      reason << "Expected value " << expected << " != " << actual;
      failed = reason.str();
    }
    else // cannot output value
      failed = "Expected value not equal";
  }
  test(file, line, failed, throw_on_fail);
}

} // namespace unit_tests
