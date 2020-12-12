#include "physicalelements.hpp"

/** @brief Constructor of FEPhysicalElement : triangle, 3 points.
 * @relates Tri3N
 */
template <>
void
Tri3N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::TRIANGLE;
    m_type_vtk      = VTK_CELL_TYPE::VTK_TRIANGLE;
    m_volref        = 0.5;

    // Points
    AddPointDefinition (0., 0., 0.);
    AddPointDefinition (1., 0., 0.);
    AddPointDefinition (0., 1., 0.);

    // Edges
    AddEdgeElement ({0, 1}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {0., -1., 0.}, {1., 0., 0.});
    AddEdgeElement ({1, 2}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {1., 1., 0.}, {-1., 1., 0.});
    AddEdgeElement ({2, 0}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {-1., 0., 0.}, {0., -1., 0.});

    return;
}
