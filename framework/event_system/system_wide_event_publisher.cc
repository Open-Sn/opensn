#include "framework/event_system/system_wide_event_publisher.h"

namespace chi
{

SystemWideEventPublisher::SystemWideEventPublisher() : chi::EventPublisher("SystemWide")
{
}

SystemWideEventPublisher&
SystemWideEventPublisher::GetInstance()
{
  static SystemWideEventPublisher instance;

  return instance;
}

void
SystemWideEventPublisher::PublishEvent(const chi::Event& event)
{
  EventPublisher::PublishEvent(event);
}

} // namespace chi
