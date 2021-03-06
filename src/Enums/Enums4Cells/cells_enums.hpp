#ifndef SRC_ENUMS_ENUMS4CELLS_CELLS_ENUMS_HPP
#define SRC_ENUMS_ENUMS4CELLS_CELLS_ENUMS_HPP

#include "../common_head.hpp"

enum class CAT_CELL_EDGE : ul_t
{
    CELL  = 0,
    EDGE  = 1,
    FIRST = CELL,
    LAST  = EDGE
};

enum class PHYSICAL_CELL_TYPE : ul_t
{
    EMPTY       = 0,
    VERTEX      = 1,
    LINE        = 2,
    TRIANGLE    = 3,
    QUADRANGLE  = 4,
    TETRAHEDRON = 5,
    FIRST       = EMPTY,
    LAST        = TETRAHEDRON
};

enum class VTK_CELL_TYPE : ul_t
{
    // Linear cells
    VTK_EMPTY_CELL       = 0,
    VTK_VERTEX           = 1,
    VTK_POLY_VERTEX      = 2,
    VTK_LINE             = 3,
    VTK_POLY_LINE        = 4,
    VTK_TRIANGLE         = 5,
    VTK_TRIANGLE_STRIP   = 6,
    VTK_POLYGON          = 7,
    VTK_PIXEL            = 8,
    VTK_QUAD             = 9,
    VTK_TETRA            = 10,
    VTK_VOXEL            = 11,
    VTK_HEXAHEDRON       = 12,
    VTK_WEDGE            = 13,
    VTK_PYRAMID          = 14,
    VTK_PENTAGONAL_PRISM = 15,
    VTK_HEXAGONAL_PRISM  = 16,

    // Quadratic, isoparametric cells
    VTK_QUADRATIC_EDGE                   = 21,
    VTK_QUADRATIC_TRIANGLE               = 22,
    VTK_QUADRATIC_QUAD                   = 23,
    VTK_QUADRATIC_POLYGON                = 36,
    VTK_QUADRATIC_TETRA                  = 24,
    VTK_QUADRATIC_HEXAHEDRON             = 25,
    VTK_QUADRATIC_WEDGE                  = 26,
    VTK_QUADRATIC_PYRAMID                = 27,
    VTK_BIQUADRATIC_QUAD                 = 28,
    VTK_TRIQUADRATIC_HEXAHEDRON          = 29,
    VTK_QUADRATIC_LINEAR_QUAD            = 30,
    VTK_QUADRATIC_LINEAR_WEDGE           = 31,
    VTK_BIQUADRATIC_QUADRATIC_WEDGE      = 32,
    VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33,
    VTK_BIQUADRATIC_TRIANGLE             = 34,

    // Cubic, isoparametric cell
    VTK_CUBIC_LINE = 35,

    // Special class of cells formed by convex group of points
    VTK_CONVEX_POINT_SET = 41,

    // Polyhedron cell (consisting of polygonal faces)
    VTK_POLYHEDRON = 42,

    // Higher order cells in parametric form
    VTK_PARAMETRIC_CURVE        = 51,
    VTK_PARAMETRIC_SURFACE      = 52,
    VTK_PARAMETRIC_TRI_SURFACE  = 53,
    VTK_PARAMETRIC_QUAD_SURFACE = 54,
    VTK_PARAMETRIC_TETRA_REGION = 55,
    VTK_PARAMETRIC_HEX_REGION   = 56,

    // Higher order cells
    VTK_HIGHER_ORDER_EDGE        = 60,
    VTK_HIGHER_ORDER_TRIANGLE    = 61,
    VTK_HIGHER_ORDER_QUAD        = 62,
    VTK_HIGHER_ORDER_POLYGON     = 63,
    VTK_HIGHER_ORDER_TETRAHEDRON = 64,
    VTK_HIGHER_ORDER_WEDGE       = 65,
    VTK_HIGHER_ORDER_PYRAMID     = 66,
    VTK_HIGHER_ORDER_HEXAHEDRON  = 67,

    // Arbitrary order Lagrange elements (formulated separated from generic higher order cells)
    VTK_LAGRANGE_CURVE         = 68,
    VTK_LAGRANGE_TRIANGLE      = 69,
    VTK_LAGRANGE_QUADRILATERAL = 70,
    VTK_LAGRANGE_TETRAHEDRON   = 71,
    VTK_LAGRANGE_HEXAHEDRON    = 72,
    VTK_LAGRANGE_WEDGE         = 73,
    VTK_LAGRANGE_PYRAMID       = 74,

    // Arbitrary order Bezier elements (formulated separated from generic higher order cells)
    VTK_BEZIER_CURVE         = 75,
    VTK_BEZIER_TRIANGLE      = 76,
    VTK_BEZIER_QUADRILATERAL = 77,
    VTK_BEZIER_TETRAHEDRON   = 78,
    VTK_BEZIER_HEXAHEDRON    = 79,
    VTK_BEZIER_WEDGE         = 80,
    VTK_BEZIER_PYRAMID       = 81,
    FIRST                    = VTK_EMPTY_CELL,
    LAST                     = VTK_BEZIER_PYRAMID
};

enum class GMSH_CELL_TYPE : ul_t
{
    GMSH_0_NODE_EMPTY                  = 0,
    GMSH_2_NODE_LINE                   = 1,
    GMSH_3_NODE_TRIANGLE               = 2,
    GMSH_4_NODE_QUADRANGLE             = 3,
    GMSH_4_NODE_TETRAHEDRON            = 4,
    GMSH_8_NODE_HEXAEDRON              = 5,
    GMSH_6_NODE_PRISM                  = 6,
    GMSH_5_NODE_PYRAMID                = 7,
    GMSH_3_NODE_QUADRATIC_LINE         = 8,
    GMSH_6_NODE_QUADRATIC_TRIANGLE     = 9,
    GMSH_9_NODE_QUADRATIC_QUADRANGLE   = 10,
    GMSH_10_NODE_QUADRATIC_TETRAHEDRON = 11,
    GMSH_27_NODE_QUADRATIC_HEXAHEDRON  = 12,
    GMSH_18_NODE_QUADRATIC_PRISM       = 13,
    GMSH_14_NODE_QUADRATIC_PYRAMID     = 14,
    GMSH_1_NODE_POINT                  = 15,
    GMSH_8_NODE_QUADRATIC_QUADRANGLE   = 16,
    GMSH_20_NODE_QUADRATIC_HEXAHEDRON  = 17,
    FIRST                              = GMSH_0_NODE_EMPTY,
    LAST                               = GMSH_20_NODE_QUADRATIC_HEXAHEDRON
};

template <>
std::string to_string (const PHYSICAL_CELL_TYPE & type);

template <>
std::string to_string (const GMSH_CELL_TYPE & type);

template <>
std::string to_string (const VTK_CELL_TYPE & type);

template <>
GMSH_CELL_TYPE Convert (const VTK_CELL_TYPE & type);

template <>
VTK_CELL_TYPE Convert (const GMSH_CELL_TYPE & type);

#endif /* SRC_ENUMS_ENUMS4CELLS_CELLS_ENUMS_HPP */
