/** @file ver1n1ddl.cpp */

#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_VERTEX, 1, 1>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_VERTEX, 1>
 */
template<>
void Ver1N1DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P0);
    m_dim = 0;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0]       = LAMBDA_EVALDOUBLE(1.);
    m_grad_phi [0]  = LAMBDA_EVALPOINT(0., 0., 0.);

    return;
}
