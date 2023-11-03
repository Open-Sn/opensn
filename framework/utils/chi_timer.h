#pragma once

#include <string>
#include <chrono>

//################################################################### CLASS DEF
namespace chi
{
/** Timer object.*/
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_;

public:
  // 00
  Timer() noexcept;
  // 01
  void Reset();
  double GetTime() const;
  std::string GetTimeString() const;
  static std::string GetLocalDateTimeString();
};

/**Puts the current thread to sleep.
 * \param time Time to sleep for.
 *
 * \note To specify different times `std::chrono` allows
 * you to change the unit with, e.g.,
 * `chi::Sleep(std::chrono::milliseconds(100))` sleeps for 100 milliseconds,
 * `std::Sleep(std::chrono::seconds(1))` sleeps for 1 second.*/
void Sleep(std::chrono::duration<double> time);
} // namespace chi


