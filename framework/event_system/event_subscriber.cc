// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/event_system/event_subscriber.h"

namespace opensn
{

void
EventSubscriber::ReceiveEventUpdate(const Event& event)
{
  // Default behavior is to not respond to events.
}

} // namespace opensn
