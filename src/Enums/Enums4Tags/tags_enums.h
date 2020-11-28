#ifndef TAGS_ENUMS_H
#define TAGS_ENUMS_H

#include "../common_head.h"

enum class INTER
{
    UNKNOWN = 0x0,
    IN   = 0x1,
    OUT   = 0x2,
    MIXED  = 0x3,
    DEFAULT = UNKNOWN,
    FIRST  = UNKNOWN,
    LAST  = MIXED
};

/*
 * NONE = 0x0,
 * WALL = 0x1,
 * INLET = 0x2,
 * WALL | INLET = 0x3,
 * OUTLET = 0x4,
 * WALL | OUTLET = 0x5,
 * INLET | OUTLET = 0x6,
 * WALL | INLET | OUTLET = 0x7,
 * DOMAIN = 0x8,
 * WALL | DOMAIN = 0x9,
 * INLET | DOMAIN = 0xa,
 * OUTLET | DOMAIN = 0xc,
 * WALL | INLET | DOMAIN = 0xb,
 * WALL | OUTLET | DOMAIN = 0xd,
 * INLET | OUTLET | DOMAIN = 0xe,
 * WALL | INLET | OUTLET | DOMAIN = 0xf,
 * FIRST = NONE,
 * LAST = DOMAIN
 */

enum class PHYS
{
    NONE  = 0x0,
    WALL  = 0x1,
    INLET  = 0x2,
    OUTLET = 0x4,
    DOMAIN = 0x8,
    DEFAULT = NONE,
    FIRST  = NONE,
    LAST  = DOMAIN
};

template<>
PHYS from_string (const std::string& s);

template<>
std::string to_string(const PHYS& type);

template<>
INTER from_string (const std::string& s);

template<>
std::string to_string(const INTER& type);

#endif // TAGS_ENUMS_H
