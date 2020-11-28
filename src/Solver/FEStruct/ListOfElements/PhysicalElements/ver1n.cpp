#include <Solver/FEStruct/ListOfElements/PhysicalElements/physicalelements.h>

/** @brief Constructor of FEPhysicalElement : vertex, 1 point.
 * @relates Ver1N
 */
template<>
void Ver1N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::VERTEX;
    m_type_vtk = VTK_CELL_TYPE::VTK_VERTEX;
    m_volref = 0.;

    // Points
    AddPointDefinition (-1., 0., 0.);

    return;
}
