#ifndef ITEM_SOLVER_TYPE_H
#define ITEM_SOLVER_TYPE_H

#include <string>
#include <ostream>

enum class ITEM_SOLVER_TYPE
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

std::string ToString (ITEM_SOLVER_TYPE type);
std::ostream& operator<< (std::ostream& out, ITEM_SOLVER_TYPE type);


#endif // ITEM_SOLVER_TYPE_H
