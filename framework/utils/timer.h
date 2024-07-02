// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <chrono>

namespace opensn
{

/** Timer object.*/
class Timer
{
private:
  std::chrono::steady_clock::time_point start_time_;

public:
  /** Default constructor.*/
  Timer() noexcept;
  /** Resets the timer to zero.*/
  void Reset();
  /** Gets the current timer value in milliseconds.*/
  double GetTime() const;
  /**
   * Obtains a string in the format of hh:mm::ss.
   */
  std::string GetTimeString() const;
  /**
   * Obtains a string in the format YYYY-MM-DD hh:mm:ss
   */
  static std::string GetLocalDateTimeString();
};

/**Puts the current thread to sleep.
 * \param time Time to sleep for.
 *
 * \note To specify different times `std::chrono` allows
 * you to change the unit with, e.g.,
 * `opensn::Sleep(std::chrono::milliseconds(100))` sleeps for 100 milliseconds,
 * `std::Sleep(std::chrono::seconds(1))` sleeps for 1 second.*/
void Sleep(std::chrono::duration<double> time);

} // namespace opensn
