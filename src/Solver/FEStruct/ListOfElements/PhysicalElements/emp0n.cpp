#include "physicalelements.hpp"

/** @brief Constructor of FEPhysicalElement : line, 2 points.
 * @relates Emp0N
 */
template <>
void
Emp0N::Build ()
{
    m_type_vtk      = VTK_CELL_TYPE::VTK_EMPTY_CELL;
    m_type_physical = PHYSICAL_CELL_TYPE::EMPTY;
    m_coor          = {};
    m_edges         = {};
    m_volref        = 0.;

    return;
}
