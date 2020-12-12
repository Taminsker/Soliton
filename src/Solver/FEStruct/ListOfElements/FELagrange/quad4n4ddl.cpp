/** @file quad4n4ddl.cpp */

#include "felagrange.hpp"

/** @brief Implemented in FELagrange<TYPE_QUADRANGLE, 4, 4>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_QUADRANGLE, 4>
 */
template <>
void
Quad4N4DDL::Build ()
{
    m_order = static_cast<int> (LAGRANGE_ORDER::ORDER_P1);
    m_dim   = 2;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE (0.25 * (1. - p->x) * (1. - p->y));
    m_phi [1] = LAMBDA_EVALDOUBLE (0.25 * (1. + p->x) * (1. - p->y));
    m_phi [2] = LAMBDA_EVALDOUBLE (0.25 * (1. + p->x) * (1. + p->y));
    m_phi [3] = LAMBDA_EVALDOUBLE (0.25 * (1. - p->x) * (1. + p->y));

    m_grad_phi [0] = LAMBDA_EVALPOINT (0.25 * (p->y - 1.), 0.25 * (p->x - 1.), 0.);
    m_grad_phi [1] = LAMBDA_EVALPOINT (0.25 * (1. - p->y), 0.25 * (-1. - p->x), 0.);
    m_grad_phi [2] = LAMBDA_EVALPOINT (0.25 * (1. + p->y), 0.25 * (1. + p->x), 0.);
    m_grad_phi [3] = LAMBDA_EVALPOINT (0.25 * (-1. - p->y), 0.25 * (1. - p->x), 0.);

    return;
}
