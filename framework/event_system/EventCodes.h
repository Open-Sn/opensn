#ifndef CHITECH_EVENTCODES_H
#define CHITECH_EVENTCODES_H

#include <string>

namespace chi
{

/**Gets the standard even code associated with the given name. If no code
* is found then 0 (i.e. Unknown Event) is returned.*/
int GetStandardEventCode(const std::string& event_name);

}

#endif // CHITECH_EVENTCODES_H
