#ifndef TAGPHYSICAL_H
#define TAGPHYSICAL_H

#include <string>

enum class TAG_PHYSICAL
{
    TAG_NONE        = 0,
    TAG_WALL        = 1,
    TAG_INLET       = 2,
    TAG_OUTLET      = 3,
    TAG_DOMAIN      = 4,
    FIRST           = TAG_WALL,
    LAST            = TAG_DOMAIN
};


std::string ToString (TAG_PHYSICAL tag);
std::ostream& operator<< (std::ostream& out, TAG_PHYSICAL tag);

#endif // TAGPHYSICAL_H
