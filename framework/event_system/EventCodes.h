#pragma once

#include <string>

namespace chi
{

/**Gets the standard even code associated with the given name. If no code
 * is found then 0 (i.e. Unknown Event) is returned.*/
int GetStandardEventCode(const std::string& event_name);

} // namespace chi


