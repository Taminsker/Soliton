#include <Solver/FEStruct/ListOfElements/PhysicalElements/physicalelements.h>

/** @brief Constructor of FEPhysicalElement : triangle, 6 points.
 * @relates Tri6N
 */
template<>
void Tri6N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::TRIANGLE;
    m_type_vtk = VTK_CELL_TYPE::VTK_QUADRATIC_TRIANGLE;
    m_volref = 0.5;

    // Points
    AddPointDefinition (0., 0., 0.);
    AddPointDefinition (1., 0., 0.);
    AddPointDefinition (0., 1., 0.);
    AddPointDefinition (0.5, 0., 0.);
    AddPointDefinition (0.5, 0.5, 0.);
    AddPointDefinition (0., 0.5, 0.);

    // Edges
    AddEdgeElement ({0, 1, 3}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {0., -1., 0.}, {1., 0., 0.});
    AddEdgeElement ({1, 2, 4}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {1., 1., 0.}, {-1., 1., 0.});
    AddEdgeElement ({2, 0, 5}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {-1., 0., 0.}, {0., -1., 0.});

    return;
}

