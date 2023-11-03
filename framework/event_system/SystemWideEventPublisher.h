#pragma once

#include "event_system/EventPublisher.h"

namespace chi
{

class SystemWideEventPublisher : public chi::EventPublisher
{
public:
  static SystemWideEventPublisher& GetInstance();

  SystemWideEventPublisher(const SystemWideEventPublisher&) = delete; // Deleted copy constructor
  SystemWideEventPublisher
  operator=(const SystemWideEventPublisher&) = delete; // Deleted assignment operator

  void PublishEvent(const chi::Event& event) override;

private:
  SystemWideEventPublisher();
};

} // namespace chi


