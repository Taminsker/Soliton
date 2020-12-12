#ifndef SRC_ENUMS_ENUMS4TAGS_TAGS_ENUMS_HPP
#define SRC_ENUMS_ENUMS4TAGS_TAGS_ENUMS_HPP

#include "../common_head.hpp"

enum class INTER : ul_t
{
    UNKNOWN = 0x0,
    IN      = 0x1,
    OUT     = 0x2,
    MIXED   = 0x3,
    DEFAULT = UNKNOWN,
    FIRST   = UNKNOWN,
    LAST    = MIXED
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

enum class PHYS : ul_t
{
    NONE    = 0x0,
    WALL    = 0x1,
    INLET   = 0x2,
    OUTLET  = 0x4,
    DOMAIN  = 0x8,
    DEFAULT = NONE,
    FIRST   = NONE,
    LAST    = DOMAIN
};

template <>
PHYS from_string (const std::string & s);

template <>
std::string to_string (const PHYS & type);

template <>
INTER from_string (const std::string & s);

template <>
std::string to_string (const INTER & type);

#endif /* SRC_ENUMS_ENUMS4TAGS_TAGS_ENUMS_HPP */
