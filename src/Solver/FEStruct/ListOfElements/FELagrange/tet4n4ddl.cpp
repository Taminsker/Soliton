/** @file tet4n4ddl.cpp */

#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_TETRAHEDRON, 4, 4>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_TETRAHEDRON, 4>
 */
template<>
void Tet4N4DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P1);
    m_dim = 2;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE(1. - p->x - p->y - p->z);
    m_phi [1] = LAMBDA_EVALDOUBLE(p->x);
    m_phi [2] = LAMBDA_EVALDOUBLE(p->y);
    m_phi [3] = LAMBDA_EVALDOUBLE(p->z);

    m_grad_phi [0] = LAMBDA_EVALPOINT(-1., -1., -1.);
    m_grad_phi [1] = LAMBDA_EVALPOINT(1., 0., 0.);
    m_grad_phi [2] = LAMBDA_EVALPOINT(0., 1., 0.);
    m_grad_phi [3] = LAMBDA_EVALPOINT(0., 0., 1.);

    return;
}
