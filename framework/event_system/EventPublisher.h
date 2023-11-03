#pragma once

#include <vector>
#include <memory>
#include <string>

namespace chi
{
class Event;
class EventSubscriber;
} // namespace chi

namespace chi
{

/**Base class for event publishers.*/
class EventPublisher
{
public:
  /**Publish the given event.*/
  virtual void PublishEvent(const chi::Event& event);
  /**Adds a subscriber to the publisher.*/
  void AddSubscriber(std::shared_ptr<chi::EventSubscriber>& subscriber_sptr);

  virtual ~EventPublisher() = default;

protected:
  explicit EventPublisher(const std::string& name);

protected:
  const std::string publisher_name_;
  std::vector<std::weak_ptr<chi::EventSubscriber>> subscribers_;
};

} // namespace chi


