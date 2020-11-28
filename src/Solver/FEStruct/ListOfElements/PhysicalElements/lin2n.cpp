#include <Solver/FEStruct/ListOfElements/PhysicalElements/physicalelements.h>

/** @brief Constructor of FEPhysicalElement : line, 2 points.
 * @relates Lin2N
 */
template<>
void Lin2N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::LINE;
    m_type_vtk = VTK_CELL_TYPE::VTK_LINE;
    m_volref = 2.;

    // Points
    AddPointDefinition (-1., 0., 0.);
    AddPointDefinition (1., 0., 0.);

    // Edges
    AddEdgeElement ({0}, GMSH_CELL_TYPE::GMSH_1_NODE_POINT, {-1., 0., 0.}, {0., 0., 1.});
    AddEdgeElement ({1}, GMSH_CELL_TYPE::GMSH_1_NODE_POINT, {1., 0., 0.}, {0., 0., -1.});

    return;
}
