/** @file lin3n3ddl.cpp */

#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_LINE, 3, 3>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_LINE, 3>
 */
template<>
void Lin3N3DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P2);
    m_dim = 1;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE(0.5 * (- p->x) * (1. - p->x));
    m_phi [1] = LAMBDA_EVALDOUBLE(0.5 * p->x * (1. + p->x));
    m_phi [2] = LAMBDA_EVALDOUBLE(0.5 * 2. * (1. - p->x * p->x));

    m_grad_phi [0] = LAMBDA_EVALPOINT(0.5 * (-1. + 2. * p->x), 0, 0);
    m_grad_phi [1] = LAMBDA_EVALPOINT(0.5 * (1. + 2. * p->x), 0, 0);
    m_grad_phi [2] = LAMBDA_EVALPOINT(0.5 * (-4) * p->x, 0, 0);

    return;
}
