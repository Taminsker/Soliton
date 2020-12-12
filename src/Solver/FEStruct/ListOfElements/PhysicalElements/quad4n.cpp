#include "physicalelements.hpp"

/** @brief Constructor of FEPhysicalElement : quadrangle, 4 point.
 * @relates Quad4N
 */
template <>
void
Quad4N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::QUADRANGLE;
    m_type_vtk      = VTK_CELL_TYPE::VTK_QUAD;
    m_volref        = 4;

    // Points
    AddPointDefinition (-1., -1., 0.);
    AddPointDefinition (1., -1., 0.);
    AddPointDefinition (1., 1., 0.);
    AddPointDefinition (-1., 1., 0.);

    // Edges
    AddEdgeElement ({0, 1}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {0., -1., 0.}, {1., 0., 0.});
    AddEdgeElement ({1, 2}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {1., 0., 0.}, {0., 1., 0.});
    AddEdgeElement ({2, 3}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {0., 1., 0.}, {-1., 0., 0.});
    AddEdgeElement ({3, 0}, GMSH_CELL_TYPE::GMSH_2_NODE_LINE, {-1., 0., 0.}, {0., -1., 0.});

    return;
}
