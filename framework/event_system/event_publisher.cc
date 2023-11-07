#include "framework/event_system/event_publisher.h"

#include "framework/event_system/event.h"
#include "framework/event_system/event_subscriber.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include <algorithm>

namespace chi
{

EventPublisher::EventPublisher(const std::string& name) : publisher_name_(name)
{
}

void
EventPublisher::PublishEvent(const chi::Event& event)
{
  size_t subs = 0;
  for (auto& subscriber_wptr : subscribers_)
    if (auto subscriber_sptr = subscriber_wptr.lock())
    {
      subscriber_sptr->ReceiveEventUpdate(event);
      ++subs;
    }
  if (Chi::log.GetVerbosity() >= 1)
    Chi::log.Log0Verbose1() << publisher_name_ << " published event name \"" << event.Name()
                            << "\"";
}

void
EventPublisher::AddSubscriber(std::shared_ptr<chi::EventSubscriber>& subscriber_sptr)
{
  std::weak_ptr<chi::EventSubscriber> wptr = subscriber_sptr;

  auto it = std::find_if(subscribers_.begin(),
                         subscribers_.end(),
                         [&wptr](const std::weak_ptr<chi::EventSubscriber>& ptr1)
                         { return ptr1.lock() == wptr.lock(); });

  if (it == subscribers_.end()) subscribers_.push_back(std::move(wptr));
}

} // namespace chi
