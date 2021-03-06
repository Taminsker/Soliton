#include "cells_enums.hpp"

#include "../def_macros.hpp"

TO_STRING_SPECIALIZATION (PHYSICAL_CELL_TYPE)
{
    CASE (VERTEX);
    CASE (LINE);
    CASE (TRIANGLE);
    CASE (QUADRANGLE);
    CASE (TETRAHEDRON);
    END_CASE_DEFAULT (EMPTY);
}

TO_STRING_SPECIALIZATION (VTK_CELL_TYPE)
{
    CASE (VTK_VERTEX);
    CASE (VTK_POLY_VERTEX);
    CASE (VTK_LINE);
    CASE (VTK_POLY_LINE);
    CASE (VTK_TRIANGLE);
    CASE (VTK_TRIANGLE_STRIP);
    CASE (VTK_POLYGON);
    CASE (VTK_PIXEL);
    CASE (VTK_QUAD);
    CASE (VTK_TETRA);
    CASE (VTK_VOXEL);
    CASE (VTK_HEXAHEDRON);
    CASE (VTK_WEDGE);
    CASE (VTK_PYRAMID);
    CASE (VTK_PENTAGONAL_PRISM);
    CASE (VTK_HEXAGONAL_PRISM);
    CASE (VTK_QUADRATIC_EDGE);
    CASE (VTK_QUADRATIC_TRIANGLE);
    CASE (VTK_QUADRATIC_QUAD);
    CASE (VTK_QUADRATIC_POLYGON);
    CASE (VTK_QUADRATIC_TETRA);
    CASE (VTK_QUADRATIC_HEXAHEDRON);
    CASE (VTK_QUADRATIC_WEDGE);
    CASE (VTK_QUADRATIC_PYRAMID);
    CASE (VTK_BIQUADRATIC_QUAD);
    CASE (VTK_TRIQUADRATIC_HEXAHEDRON);
    CASE (VTK_QUADRATIC_LINEAR_QUAD);
    CASE (VTK_QUADRATIC_LINEAR_WEDGE);
    CASE (VTK_BIQUADRATIC_QUADRATIC_WEDGE);
    CASE (VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON);
    CASE (VTK_BIQUADRATIC_TRIANGLE);
    CASE (VTK_CUBIC_LINE);
    CASE (VTK_CONVEX_POINT_SET);
    CASE (VTK_POLYHEDRON);
    CASE (VTK_PARAMETRIC_CURVE);
    CASE (VTK_PARAMETRIC_SURFACE);
    CASE (VTK_PARAMETRIC_TRI_SURFACE);
    CASE (VTK_PARAMETRIC_QUAD_SURFACE);
    CASE (VTK_PARAMETRIC_TETRA_REGION);
    CASE (VTK_PARAMETRIC_HEX_REGION);
    CASE (VTK_HIGHER_ORDER_EDGE);
    CASE (VTK_HIGHER_ORDER_TRIANGLE);
    CASE (VTK_HIGHER_ORDER_QUAD);
    CASE (VTK_HIGHER_ORDER_POLYGON);
    CASE (VTK_HIGHER_ORDER_TETRAHEDRON);
    CASE (VTK_HIGHER_ORDER_WEDGE);
    CASE (VTK_HIGHER_ORDER_PYRAMID);
    CASE (VTK_HIGHER_ORDER_HEXAHEDRON);
    CASE (VTK_LAGRANGE_CURVE);
    CASE (VTK_LAGRANGE_TRIANGLE);
    CASE (VTK_LAGRANGE_QUADRILATERAL);
    CASE (VTK_LAGRANGE_TETRAHEDRON);
    CASE (VTK_LAGRANGE_HEXAHEDRON);
    CASE (VTK_LAGRANGE_WEDGE);
    CASE (VTK_LAGRANGE_PYRAMID);
    CASE (VTK_BEZIER_CURVE);
    CASE (VTK_BEZIER_TRIANGLE);
    CASE (VTK_BEZIER_QUADRILATERAL);
    CASE (VTK_BEZIER_TETRAHEDRON);
    CASE (VTK_BEZIER_HEXAHEDRON);
    CASE (VTK_BEZIER_WEDGE);
    CASE (VTK_BEZIER_PYRAMID);
    END_CASE_DEFAULT (VTK_EMPTY_CELL);
}

TO_STRING_SPECIALIZATION (GMSH_CELL_TYPE)
{
    CASE (GMSH_2_NODE_LINE);
    CASE (GMSH_3_NODE_TRIANGLE);
    CASE (GMSH_4_NODE_QUADRANGLE);
    CASE (GMSH_4_NODE_TETRAHEDRON);
    CASE (GMSH_8_NODE_HEXAEDRON);
    CASE (GMSH_6_NODE_PRISM);
    CASE (GMSH_5_NODE_PYRAMID);
    CASE (GMSH_3_NODE_QUADRATIC_LINE);
    CASE (GMSH_6_NODE_QUADRATIC_TRIANGLE);
    CASE (GMSH_9_NODE_QUADRATIC_QUADRANGLE);
    CASE (GMSH_10_NODE_QUADRATIC_TETRAHEDRON);
    CASE (GMSH_27_NODE_QUADRATIC_HEXAHEDRON);
    CASE (GMSH_18_NODE_QUADRATIC_PRISM);
    CASE (GMSH_14_NODE_QUADRATIC_PYRAMID);
    CASE (GMSH_1_NODE_POINT);
    CASE (GMSH_8_NODE_QUADRATIC_QUADRANGLE);
    CASE (GMSH_20_NODE_QUADRATIC_HEXAHEDRON);
    END_CASE_DEFAULT (GMSH_0_NODE_EMPTY);
}

CONVERT_SPECIALIZATION (VTK_CELL_TYPE, GMSH_CELL_TYPE)
{
    CASE_CONVERT (VTK_LINE, GMSH_2_NODE_LINE);
    CASE_CONVERT (VTK_TRIANGLE, GMSH_3_NODE_TRIANGLE);
    CASE_CONVERT (VTK_QUAD, GMSH_4_NODE_QUADRANGLE);
    CASE_CONVERT (VTK_TETRA, GMSH_4_NODE_TETRAHEDRON);
    CASE_CONVERT (VTK_HEXAHEDRON, GMSH_8_NODE_HEXAEDRON);
    //  CASE_CONVERT(VTK_PYRAMID,          GMSH_6_NODE_PRISM);
    CASE_CONVERT (VTK_QUADRATIC_PYRAMID, GMSH_5_NODE_PYRAMID);
    CASE_CONVERT (VTK_QUADRATIC_EDGE, GMSH_3_NODE_QUADRATIC_LINE);
    CASE_CONVERT (VTK_QUADRATIC_TRIANGLE, GMSH_6_NODE_QUADRATIC_TRIANGLE);
    CASE_CONVERT (VTK_QUADRATIC_QUAD, GMSH_9_NODE_QUADRATIC_QUADRANGLE);
    CASE_CONVERT (VTK_QUADRATIC_TETRA, GMSH_10_NODE_QUADRATIC_TETRAHEDRON);
    CASE_CONVERT (VTK_QUADRATIC_HEXAHEDRON, GMSH_27_NODE_QUADRATIC_HEXAHEDRON);
    //  CASE_CONVERT(VTK_EMPTY_CELL,        GMSH_18_NODE_QUADRATIC_PRISM);
    //  CASE_CONVERT(VTK_QUADRATIC_PYRAMID,     GMSH_14_NODE_QUADRATIC_PYRAMID);
    CASE_CONVERT (VTK_VERTEX, GMSH_1_NODE_POINT);
    CASE_CONVERT (VTK_BIQUADRATIC_QUAD, GMSH_8_NODE_QUADRATIC_QUADRANGLE);
    CASE_CONVERT (VTK_TRIQUADRATIC_HEXAHEDRON, GMSH_20_NODE_QUADRATIC_HEXAHEDRON);
    END_CASE_DEFAULT_CONVERT (GMSH_0_NODE_EMPTY);
}

CONVERT_SPECIALIZATION (GMSH_CELL_TYPE, VTK_CELL_TYPE)
{
    CASE_CONVERT (GMSH_2_NODE_LINE, VTK_LINE);
    CASE_CONVERT (GMSH_3_NODE_TRIANGLE, VTK_TRIANGLE);
    CASE_CONVERT (GMSH_4_NODE_QUADRANGLE, VTK_QUAD);
    CASE_CONVERT (GMSH_4_NODE_TETRAHEDRON, VTK_TETRA);
    CASE_CONVERT (GMSH_8_NODE_HEXAEDRON, VTK_HEXAHEDRON);
    //  CASE_CONVERT(GMSH_6_NODE_PRISM,           VTK_PYRAMID);
    CASE_CONVERT (GMSH_5_NODE_PYRAMID, VTK_QUADRATIC_PYRAMID);
    CASE_CONVERT (GMSH_3_NODE_QUADRATIC_LINE, VTK_QUADRATIC_EDGE);
    CASE_CONVERT (GMSH_6_NODE_QUADRATIC_TRIANGLE, VTK_QUADRATIC_TRIANGLE);
    CASE_CONVERT (GMSH_9_NODE_QUADRATIC_QUADRANGLE, VTK_QUADRATIC_QUAD);
    CASE_CONVERT (GMSH_10_NODE_QUADRATIC_TETRAHEDRON, VTK_QUADRATIC_TETRA);
    CASE_CONVERT (GMSH_27_NODE_QUADRATIC_HEXAHEDRON, VTK_QUADRATIC_HEXAHEDRON);
    //  CASE_CONVERT(GMSH_18_NODE_QUADRATIC_PRISM,     VTK_EMPTY_CELL);
    CASE_CONVERT (GMSH_14_NODE_QUADRATIC_PYRAMID, VTK_QUADRATIC_PYRAMID);
    CASE_CONVERT (GMSH_1_NODE_POINT, VTK_VERTEX);
    CASE_CONVERT (GMSH_8_NODE_QUADRATIC_QUADRANGLE, VTK_BIQUADRATIC_QUAD);
    CASE_CONVERT (GMSH_20_NODE_QUADRATIC_HEXAHEDRON, VTK_TRIQUADRATIC_HEXAHEDRON);
    END_CASE_DEFAULT_CONVERT (VTK_EMPTY_CELL);
}

#include "../undef_macros.hpp"
