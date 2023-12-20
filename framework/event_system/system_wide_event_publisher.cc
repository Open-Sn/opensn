#include "framework/event_system/system_wide_event_publisher.h"

namespace opensn
{

SystemWideEventPublisher::SystemWideEventPublisher() : EventPublisher("SystemWide")
{
}

SystemWideEventPublisher&
SystemWideEventPublisher::GetInstance()
{
  static SystemWideEventPublisher instance;

  return instance;
}

void
SystemWideEventPublisher::PublishEvent(const Event& event)
{
  EventPublisher::PublishEvent(event);
}

} // namespace opensn
