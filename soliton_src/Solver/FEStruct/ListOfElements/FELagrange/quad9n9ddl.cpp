/** @file quad9n9ddl.cpp */

#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_QUADRANGLE, 9, 9>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_QUADRANGLE, 9>
 */
template<>
void Quad9N9DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P2);
    m_dim = 2;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE(0.25 * p->x * p->y * (1. - p->x) * (1 - p->y));
    m_phi [1] = LAMBDA_EVALDOUBLE(- 0.25 * p->x * p->y * (1. + p->x) * (1 - p->y));
    m_phi [2] = LAMBDA_EVALDOUBLE(0.25 * p->x * p->y * (1. + p->x) * (1 + p->y));
    m_phi [3] = LAMBDA_EVALDOUBLE(- 0.25 * p->x * p->y * (1. - p->x) * (1 + p->y));
    m_phi [4] = LAMBDA_EVALDOUBLE(- 0.5 * p->y * (1. - p->x * p->x) * (1. - p->y));
    m_phi [5] = LAMBDA_EVALDOUBLE(0.5 * p->x * (1. + p->x) * (1. - p->y * p->y));
    m_phi [6] = LAMBDA_EVALDOUBLE(0.5 * p->y * (1. - p->x * p->x) * (1. + p->y));
    m_phi [7] = LAMBDA_EVALDOUBLE(- 0.5 * p->x * (1. - p->x) * (1. - p->y * p->y));
    m_phi [8] = LAMBDA_EVALDOUBLE((1. - p->x * p->x) * (1. - p->y * p->y));


    m_grad_phi [0] = LAMBDA_EVALPOINT(0.25 * p->y * (1. - 2. * p->x) * (1. - p->y),
                                   0.25 * p->x * (1. - p->x) * (1. - 2. * p->y),
                                   0.);
    m_grad_phi [1] = LAMBDA_EVALPOINT(- 0.25 * p->y * (1. + 2. * p->x) * (1. - p->y),
                                   - 0.25 * p->x * (1. + p->x) * (1. - 2. * p->y),
                                   0.);
    m_grad_phi [2] = LAMBDA_EVALPOINT(0.25 * p->y * (1. + 2. * p->x) * (1. + p->y),
                                   0.25 * p->x * (1. + p->x) * (1. + 2. * p->y),
                                   0.);
    m_grad_phi [3] = LAMBDA_EVALPOINT(- 0.25 * p->y * (1. - 2. * p->x) * (1. + p->y),
                                   - 0.25 * p->x * (1. - p->x) * (1. + 2. * p->y),
                                   0.);
    m_grad_phi [4] = LAMBDA_EVALPOINT(- 0.5 * (-2.) * (1. - p->x) * p->x * p->y,
                                   - 0.5 * (1. - p->x * p->x) * (1. - 2. * p->y),
                                   0.);
    m_grad_phi [5] = LAMBDA_EVALPOINT(0.5 * (1. + 2. * p->x) * (1. - p->y * p->y),
                                   0.5 * (-2.) * (1. + p->x) * p->x * p->y,
                                   0.);
    m_grad_phi [6] = LAMBDA_EVALPOINT(0.5 * (-2.) * (1. + p->y) * p->x * p->y,
                                   0.5 * (1. - p->x * p->x) * (1. + 2. * p->y),
                                   0.);
    m_grad_phi [7] = LAMBDA_EVALPOINT(- 0.5 * (1. - 2. * p->x) * (1. - p->y * p->y) ,
                                   - 0.5 * (-2.) * (1. - p->x) * p->x * p->y,
                                   0.);
    m_grad_phi [8] = LAMBDA_EVALPOINT(- 2 * (1. - p->y * p->y) * p->x,
                                   - 2 * (1. - p->x * p->x) * p->y,
                                   0.);
    return;
}
