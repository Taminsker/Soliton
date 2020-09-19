#include "tagphysical.h"

std::string ToString (TAG_PHYSICAL tag)
{
    switch (tag)
    {
    case TAG_PHYSICAL::TAG_NONE:
        return "TAG_NONE";
    case TAG_PHYSICAL::TAG_WALL:
        return "TAG_WALL";
    case TAG_PHYSICAL::TAG_INLET:
        return "TAG_INLET";
    case TAG_PHYSICAL::TAG_OUTLET:
        return "TAG_OUTLET";
    case TAG_PHYSICAL::TAG_DOMAIN:
        return "TAG_DOMAIN";
    }

    return "TAG_NONE";
}

std::ostream& operator<< (std::ostream& out, TAG_PHYSICAL tag)
{
    out << ToString (tag);
    return out;
}
