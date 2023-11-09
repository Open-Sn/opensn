#pragma once

#include "framework/event_system/event.h"
#include <string>

namespace chi
{

/**Gets the standard even code associated with the given name. If no code
 * is found then 0 (i.e. Unknown Event) is returned.*/
Event::EventCode GetStandardEventCode(const std::string& event_name);

} // namespace chi
