#pragma once

#include "framework/event_system/event_publisher.h"

namespace chi
{

class SystemWideEventPublisher : public chi::EventPublisher
{
public:
  static SystemWideEventPublisher& GetInstance();

  /// Deleted copy constructor
  SystemWideEventPublisher(const SystemWideEventPublisher&) = delete;
  SystemWideEventPublisher
  /// Deleted assignment operator
  operator=(const SystemWideEventPublisher&) = delete;

  void PublishEvent(const chi::Event& event) override;

private:
  SystemWideEventPublisher();
};

} // namespace chi
