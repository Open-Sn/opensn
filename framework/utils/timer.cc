#include "framework/utils/timer.h"

#include <cmath>
#include <ctime>
#include <thread>

namespace opensn
{

Timer::Timer() noexcept
{
  start_time_ = std::chrono::steady_clock::now();
}

void
Timer::Reset()
{
  start_time_ = std::chrono::steady_clock::now();
}

double
Timer::GetTime() const
{
  using namespace std::chrono;

  steady_clock::time_point newTime = std::chrono::steady_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(newTime - start_time_);

  return time_span.count() * 1000.0;
}

std::string
Timer::GetTimeString() const
{
  double time_sec = this->GetTime() / 1000.0;
  int hours = std::floor(time_sec / 60 / 60);
  int minutes = std::floor((time_sec - 60 * 60 * hours) / 60);
  int seconds = (int)time_sec - 3600 * hours - 60 * minutes;

  char buff[100];
  snprintf(buff, 100, "%02d:%02d:%02d", hours, minutes, seconds);

  return {buff};
}

std::string
Timer::GetLocalDateTimeString()
{
  using namespace std::chrono;
  std::time_t now = system_clock::to_time_t(system_clock::now());

  char s[30];
  size_t end = std::strftime(s, 30, "%Y-%m-%d %H:%M:%S", std::localtime(&now));
  s[29] = '\0';
  if (end < 30)
    s[end] = '\0';
  return s;
}

void
Sleep(std::chrono::duration<double> time)
{
  std::this_thread::sleep_for(time);
}

} // namespace opensn
