// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/event_system/event_publisher.h"

namespace opensn
{

class SystemWideEventPublisher : public EventPublisher
{
public:
  static SystemWideEventPublisher& Instance();

  /// Deleted copy constructor
  SystemWideEventPublisher(const SystemWideEventPublisher&) = delete;

  /// Deleted assignment operator
  SystemWideEventPublisher operator=(const SystemWideEventPublisher&) = delete;

  void PublishEvent(const Event& event) override;

private:
  SystemWideEventPublisher();
};

} // namespace opensn
