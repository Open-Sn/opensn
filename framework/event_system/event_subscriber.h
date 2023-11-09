#pragma once

namespace chi
{

class Event;

class EventSubscriber
{
public:
  EventSubscriber() = default;

  /**A method called by publishers to inform the object of events, only if
   * the object subscribed to the publisher.*/
  virtual void ReceiveEventUpdate(const Event& event);

  virtual ~EventSubscriber() = default;
};

} // namespace chi
