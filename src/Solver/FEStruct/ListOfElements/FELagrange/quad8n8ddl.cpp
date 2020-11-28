/** @file quad8n8ddl.cpp */

#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>

/** @brief Implemented in FELagrange<TYPE_QUADRANGLE, 8, 8>.
 * @relates FELagrange
 * @implements FEPhysicalElement<TYPE_QUADRANGLE, 8>
 */
template<>
void Quad8N8DDL::Build ()
{
    m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P2);
    m_dim = 2;

    m_phi.resize (this->m_coor.size ());
    m_grad_phi.resize (this->m_coor.size ());

    m_phi [0] = LAMBDA_EVALDOUBLE(- 0.25 * (1. - p->x) * (1. - p->y) * (1 + p->x + p->y));
    m_phi [1] = LAMBDA_EVALDOUBLE(- 0.25 * (1. + p->x) * (1. - p->y) * (1 - p->x + p->y));
    m_phi [2] = LAMBDA_EVALDOUBLE(- 0.25 * (1. + p->x) * (1. + p->y) * (1 - p->x - p->y));
    m_phi [3] = LAMBDA_EVALDOUBLE(- 0.25 * (1. - p->x) * (1. + p->y) * (1 + p->x - p->y));
    m_phi [4] = LAMBDA_EVALDOUBLE(0.5 * (1. - p->x * p->x) * (1. - p->y));
    m_phi [5] = LAMBDA_EVALDOUBLE(0.5 * (1. + p->x) * (1. - p->y * p->y));
    m_phi [6] = LAMBDA_EVALDOUBLE(0.5 * (1. - p->x * p->x) * (1. + p->y));
    m_phi [7] = LAMBDA_EVALDOUBLE(0.5 * (1. - p->x) * (1. - p->y * p->y));

    m_grad_phi [0] = LAMBDA_EVALPOINT(0.25 * (1. - p->y) * (2. * p->x + p->y),
                                      0.25 * (1. - p->x) * (p->x + 2. * p->y),
                                      0.);
    m_grad_phi [1] = LAMBDA_EVALPOINT(0.25 * (1. - p->y) * (2. * p->x - p->y),
                                      -0.25 * (1. + p->x) * (p->x - 2. * p->y),
                                      0.);
    m_grad_phi [2] = LAMBDA_EVALPOINT(0.25 * (1. + p->y) * (2. * p->x + p->y),
                                      0.25 * (1. + p->x) * (p->x + 2. * p->y),
                                      0.);
    m_grad_phi [3] = LAMBDA_EVALPOINT(0.25 * (1. + p->y) * (2. * p->x - p->y),
                                      -0.25 * (1. - p->x) * (p->x - 2. * p->y),
                                      0.);
    m_grad_phi [4] = LAMBDA_EVALPOINT(-0.5 * 2. * p->x * (1. - p->y),
                                      -0.5 * (1. - p->x * p->x),
                                      0.);
    m_grad_phi [5] = LAMBDA_EVALPOINT(0.5 * (1. - p->y * p->y),
                                      -0.5 * p->y * (1. + p->x),
                                      0.);
    m_grad_phi [6] = LAMBDA_EVALPOINT(-0.5 * p->x * (1. + p->y),
                                      0.5 * (1. - p->x * p->x),
                                      0.);
    m_grad_phi [7] = LAMBDA_EVALPOINT(-0.5 * (1. - p->y * p->y),
                                      -0.5 * p->y * (1. - p->x),
                                      0.);
    return;
}



















