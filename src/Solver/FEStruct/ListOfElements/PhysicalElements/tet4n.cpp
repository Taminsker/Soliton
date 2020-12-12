#include "physicalelements.hpp"

/** @brief Constructor of FEPhysicalElement : tetrahedron, 4 point.
 * @relates Tet4N
 */
template <>
void
Tet4N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::TETRAHEDRON;
    m_type_vtk      = VTK_CELL_TYPE::VTK_TETRA;
    m_volref        = 1. / 3.;

    // Points
    AddPointDefinition (0., 0., 0.);
    AddPointDefinition (1., 0., 0.);
    AddPointDefinition (0., 1., 0.);
    AddPointDefinition (0., 0., 1.);

    // Edges
    AddEdgeElement ({0, 1, 3}, GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE, {0., -1., 0.}, {1., 0., 0.});
    AddEdgeElement ({0, 2, 1}, GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE, {0., 0., -1.}, {0., 1., 0.});
    AddEdgeElement ({0, 3, 2}, GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE, {-1., 0., 0.}, {0., 0., 1.});
    AddEdgeElement ({1, 2, 3}, GMSH_CELL_TYPE::GMSH_3_NODE_TRIANGLE, {1., 1., 0.}, {-1., 1., 0.});

    return;
}
