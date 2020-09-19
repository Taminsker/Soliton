/** @file tri6n6ddl.cpp */
#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_TRIANGLE, 6, 6>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_TRIANGLE, 6>
 */
template <>
void Tri6N6DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P2);
    m_dim = 2;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_FUNCTIONDOUBLE(
                double lamb = 1. - p->x - p->y;
            return lamb * (2. * lamb - 1.););
    m_phi [1] = LAMBDA_EVALDOUBLE(p->x * (2. * p->x - 1.));
    m_phi [2] = LAMBDA_EVALDOUBLE(p->y * (2. * p->y - 1.));
    m_phi [3] = LAMBDA_FUNCTIONDOUBLE(
                double lamb = 1. - p->x - p->y;
            return 4. * p->x * lamb;);
    m_phi [4] = LAMBDA_EVALDOUBLE(4 * p->x * p->y);
    m_phi [5] = LAMBDA_FUNCTIONDOUBLE(
                double lamb = 1. - p->x - p->y;
            return 4. * p->y * lamb;);

    m_grad_phi [0] = LAMBDA_FUNCTIONPOINT(
                double lamb = 1. - p->x - p->y;
            return Point(1. - 4. * lamb, 1. - 4. * lamb, 0););
    m_grad_phi [1] = LAMBDA_EVALPOINT(4. * p->x - 1., 0., 0.);
    m_grad_phi [2] = LAMBDA_EVALPOINT(0., 4. * p->y - 1., 0.);

    m_grad_phi [3] = LAMBDA_FUNCTIONPOINT(
                double lamb = 1. - p->x - p->y;
            return Point(4. * (lamb - p->x), -4 * p->x, 0););

    m_grad_phi [4] = LAMBDA_EVALPOINT(4 * p->y, 4 * p->x, 0.);
    m_grad_phi [5] = LAMBDA_FUNCTIONPOINT(
                double lamb = 1. - p->x - p->y;
            return Point(-4 * p->y, 4. * (lamb - p->y), 0););

    return;
}
