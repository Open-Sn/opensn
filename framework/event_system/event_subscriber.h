// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

class Event;

class EventSubscriber
{
public:
  EventSubscriber() = default;

  /**
   * A method called by publishers to inform the object of events, only if the object subscribed to
   * the publisher.
   */
  virtual void ReceiveEventUpdate(const Event& event);

  virtual ~EventSubscriber() = default;
};

} // namespace opensn
