// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>
#include <memory>
#include <string>

namespace opensn
{
class Event;
class EventSubscriber;

/**Base class for event publishers.*/
class EventPublisher
{
public:
  /**Publish the given event.*/
  virtual void PublishEvent(const Event& event);
  /**Adds a subscriber to the publisher.*/
  void AddSubscriber(std::shared_ptr<EventSubscriber>& subscriber_sptr);

  virtual ~EventPublisher() = default;

protected:
  explicit EventPublisher(const std::string& name);

protected:
  const std::string publisher_name_;
  std::vector<std::weak_ptr<EventSubscriber>> subscribers_;
};

} // namespace opensn
