#ifndef SUPERSOLVER_ENUMS_H
#define SUPERSOLVER_ENUMS_H

#include "../common_head.h"

enum class SCH_T
{
    NO_TIME = 0x0,
    EULER_IMPLICIT,
    EULER_EXPLICIT,
    TVD_RK_2,
    TVD_RK_4,
    NEWMARK_SCHEME,
    FIRST = NO_TIME,
    LAST = NEWMARK_SCHEME,
    DEFAULT = NO_TIME
};

enum class ITEM_T
{
    /** Empty item */
    EMPTY = 0x0,
    /** <phi_i | phi_j> item */
    PHI_PHI,
    /** <chi phi_i | phi_j> item */
    DAMPING_PHI_PHI,
    /** <Nabla phi_i | Nabla phi_j > item */
    GRAD_GRAD,
    /** <f | phi_j> item */
    SECOND_MEMBER,
    /** Dirichlet Boundary with penalization item */
    DIRICHLET_BOUNDARY,
    /** Neumann Boundary with penalization item */
    NEUMANN_BOUNDARY,
    /** Dirichlet Boundary with penalization and SBM item */
    DIRICHLET_BOUNDARY_SBM,
    /** Neumann Boundary with penalization and SBM item */
    NEUMANN_BOUNDARY_SBM,
    FIRST = EMPTY,
    LAST = NEUMANN_BOUNDARY_SBM
};


template<>
SCH_T from_string (const std::string& s);

template<>
std::string to_string(const SCH_T& type);

template<>
ITEM_T from_string (const std::string& s);

template<>
std::string to_string(const ITEM_T& type);
#endif // SUPERSOLVER_ENUMS_H
