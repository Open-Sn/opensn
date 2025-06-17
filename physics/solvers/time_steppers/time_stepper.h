// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

namespace opensn
{

enum class TimeStepStatus
{
  SUCCESS = 0,
  FAILURE = 1,
  NEUTRAL = 2
};

/// Base class for all timestep controllers.
class TimeStepper : public Object
{
public:
  /// Overridable method to get the timestep size.
  double GetTimeStepSize() const;

  /// Returns the current controller time.
  double GetTime() const;

  /// Returns the current time index.
  size_t GetTimeStepIndex() const;

  /// Returns the current controller start_time.
  double GetStartTime() const;

  /// Returns the current controller end_time.
  double GetEndTime() const;

  /// Returns the current controller max time steps.
  double GetMaxTimeSteps() const;

  /// If start_time <= time <= end_time, this will return true.
  bool IsActive() const;

  /// Manually set the time step size.
  void SetTimeStepSize(double dt);

  /// Manually set the current time.
  void SetTime(double time);

  /// Manually set the start_time.
  void SetStartTime(double time);

  /// Manually set the end_time.
  void SetEndTime(double time);

  /// Manually set the maximum number of time steps. A negative number disables this check.
  void SetMaxTimeSteps(int n);

  /// Manually sets the minimum time step size.
  void SetMinimumTimeStepSize(double dt_min);

  /**
   * Advances the controller's state. The most basic action here is to advance time and the time
   * index. If the solver is at or beyond its end time then it will return false. Otherwise it will
   * advance and return true.
   */
  virtual void Advance();

  /**
   * Adapts according to the timestep status. If it could provide a change it returns true,
   * otherwise false.
   */
  virtual bool Adapt(TimeStepStatus time_step_status) { return false; }

  /// Builds a formatted string of the time information.
  std::string StringTimeInfo(bool old_time = false) const;

protected:
  static InputParameters GetInputParameters();
  explicit TimeStepper(const InputParameters& params);

  double dt_;
  double time_;
  size_t t_index_;

  double start_time_;
  double end_time_;
  int max_time_steps_;
  double dt_min_;

  const double general_tolerance_;
  /// Last dt before finishing.
  double last_dt_;
};

} // namespace opensn
