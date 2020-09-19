/** @file lin2n2ddl.cpp */

#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_LINE, 2, 2>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_LINE, 2>
 */
template<>
void Lin2N2DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P1);
    m_dim = 1;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE(0.5 * (1. - p->x));
    m_phi [1] = LAMBDA_EVALDOUBLE(0.5 * (1. + p->x));

    m_grad_phi [0] = LAMBDA_EVALPOINT(-0.5, 0, 0);
    m_grad_phi [1] = LAMBDA_EVALPOINT(0.5, 0, 0);

    return;
}
