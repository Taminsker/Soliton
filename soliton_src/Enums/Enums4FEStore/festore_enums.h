#ifndef FESTORE_ENUMS_H
#define FESTORE_ENUMS_H

#include "../common_head.h"

enum class LAGRANGE_ORDER
{
    ORDER_P0 = 0,
    ORDER_P1 = 1,
    ORDER_P2 = 2,
    ORDER_P3 = 3,
    ORDER_P4 = 4,
    FIRST   = ORDER_P0,
    LAST    = ORDER_P4
};

enum class FE_CLASS_TYPE
{
    EMPTY = 0,
    LAGRANGE = 1,
    HERMITE = 2,
    FIRST = EMPTY,
    LAST = HERMITE
};


#endif // FESTORE_ENUMS_H
