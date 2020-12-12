#include "physicalelements.hpp"

/** @brief Constructor of FEPhysicalElement : quadrangle, 8 point.
 * @relates Quad8N
 */
template <>
void
Quad8N::Build ()
{
    m_type_physical = PHYSICAL_CELL_TYPE::QUADRANGLE;
    m_type_vtk      = VTK_CELL_TYPE::VTK_BIQUADRATIC_QUAD;
    m_volref        = 1;

    // Points
    AddPointDefinition (-1., -1., 0.);
    AddPointDefinition (1., -1., 0.);
    AddPointDefinition (1., 1., 0.);
    AddPointDefinition (-1., 1., 0.);
    AddPointDefinition (0., -1., 0.);
    AddPointDefinition (1., 0., 0.);
    AddPointDefinition (0., 1., 0.);
    AddPointDefinition (-1., 0., 0.);

    // Edges
    AddEdgeElement ({0, 1, 4}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {0., -1., 0.}, {1., 0., 0.});
    AddEdgeElement ({1, 2, 5}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {1., 0., 0.}, {0., 1., 0.});
    AddEdgeElement ({2, 3, 6}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {0., 1., 0.}, {-1., 0., 0.});
    AddEdgeElement ({3, 0, 7}, GMSH_CELL_TYPE::GMSH_3_NODE_QUADRATIC_LINE, {-1., 0., 0.}, {0., -1., 0.});

    return;
}
