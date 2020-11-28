/** @file emp0n0ddl.cpp */
#include <Solver/FEStruct/ListOfElements/FELagrange/felagrange.h>


/** @brief Implemented in FELagrange<TYPE_EMPTY, 0, 0>.
 * @relates FELagrange
 */
template <>
void Emp0N0DDL::Build ()
{
    this->m_order = static_cast<int>(LAGRANGE_ORDER::ORDER_P0);
    this->m_dim = 0;
    this->m_type_physical = PHYSICAL_CELL_TYPE::EMPTY;

    return;
}
