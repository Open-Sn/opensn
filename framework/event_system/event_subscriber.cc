#include "framework/event_system/event_subscriber.h"

namespace opensn
{

void
EventSubscriber::ReceiveEventUpdate(const Event& event)
{
  // Default behavior is to not respond to events.
}

} // namespace opensn
