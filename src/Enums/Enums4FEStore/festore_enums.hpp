#ifndef SRC_ENUMS_ENUMS4FESTORE_FESTORE_ENUMS_HPP
#define SRC_ENUMS_ENUMS4FESTORE_FESTORE_ENUMS_HPP

#include "../common_head.hpp"

enum class LAGRANGE_ORDER : ul_t
{
    ORDER_P0 = 0,
    ORDER_P1 = 1,
    ORDER_P2 = 2,
    ORDER_P3 = 3,
    ORDER_P4 = 4,
    FIRST    = ORDER_P0,
    LAST     = ORDER_P4
};

enum class FE_CLASS_TYPE : ul_t
{
    EMPTY    = 0,
    LAGRANGE = 1,
    HERMITE  = 2,
    FIRST    = EMPTY,
    LAST     = HERMITE
};

#endif /* SRC_ENUMS_ENUMS4FESTORE_FESTORE_ENUMS_HPP */
