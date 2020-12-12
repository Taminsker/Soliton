/** @file tri3n3ddl.cpp */

#include "felagrange.hpp"

/** @brief Implemented in FELagrange<TYPE_TRIANGLE, 3, 3>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_TRIANGLE, 3>
 */
template <>
void
Tri3N3DDL::Build ()
{
    m_order = static_cast<int> (LAGRANGE_ORDER::ORDER_P1);
    m_dim   = 2;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE (1. - p->x - p->y);
    m_phi [1] = LAMBDA_EVALDOUBLE (p->x);
    m_phi [2] = LAMBDA_EVALDOUBLE (p->y);

    m_grad_phi [0] = LAMBDA_EVALPOINT (-1., -1., 0.);
    m_grad_phi [1] = LAMBDA_EVALPOINT (1., 0., 0.);
    m_grad_phi [2] = LAMBDA_EVALPOINT (0., 1., 0.);

    return;
}
