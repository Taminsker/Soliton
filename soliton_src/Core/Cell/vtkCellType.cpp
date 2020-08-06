#include "vtkCellType.h"


std::string GetNameVTKType (VTKCellType idtype)
{
    std::string type;
    switch (idtype) {
    case VTK_EMPTY_CELL:
        type = "VTK_EMPTY_CELL";
        break;
    case VTK_VERTEX:
        type = "VTK_VERTEX";
        break;
    case VTK_POLY_VERTEX:
        type = "VTK_POLY_VERTEX";
        break;
    case VTK_LINE:
        type = "VTK_LINE";
        break;
    case VTK_POLY_LINE:
        type = "VTK_POLY_LINE ";
        break;
    case VTK_TRIANGLE:
        type = "VTK_TRIANGLE";
        break;
    case VTK_TRIANGLE_STRIP:
        type = "VTK_TRIANGLE_STRIP";
        break;
    case VTK_POLYGON:
        type = "VTK_POLYGON";
        break;
    case VTK_PIXEL:
        type = "VTK_PIXEL";
        break;
    case VTK_QUAD:
        type = "VTK_QUAD ";
        break;
    case VTK_TETRA:
        type = "VTK_TETRA  ";
        break;
    case VTK_VOXEL:
        type = "VTK_VOXEL";
        break;
    case VTK_HEXAHEDRON:
        type = "VTK_HEXAHEDRON";
        break;
    case VTK_WEDGE:
        type = "VTK_WEDGE";
        break;
    case VTK_PYRAMID:
        type = "VTK_PYRAMID ";
        break;
    case VTK_PENTAGONAL_PRISM:
        type = "VTK_PENTAGONAL_PRISM";
        break;
    case VTK_HEXAGONAL_PRISM:
        type = "VTK_HEXAGONAL_PRISM";
        break;
    case VTK_QUADRATIC_EDGE:
        type = "VTK_QUADRATIC_EDGE";
        break;
    case VTK_QUADRATIC_TRIANGLE:
        type = "VTK_QUADRATIC_TRIANGLE ";
        break;
    case VTK_QUADRATIC_QUAD:
        type = "VTK_QUADRATIC_QUAD";
        break;
    case VTK_QUADRATIC_POLYGON:
        type = "VTK_QUADRATIC_POLYGON ";
        break;
    case VTK_QUADRATIC_TETRA:
        type = "VTK_QUADRATIC_TETRA ";
        break;
    case VTK_QUADRATIC_HEXAHEDRON:
        type = "VTK_QUADRATIC_HEXAHEDRON  ";
        break;
    case VTK_QUADRATIC_WEDGE:
        type = "VTK_QUADRATIC_WEDGE      ";
        break;
    case VTK_QUADRATIC_PYRAMID:
        type = "VTK_QUADRATIC_PYRAMID";
        break;
    case VTK_BIQUADRATIC_QUAD:
        type = "VTK_BIQUADRATIC_QUAD ";
        break;
    case VTK_TRIQUADRATIC_HEXAHEDRON:
        type = "VTK_TRIQUADRATIC_HEXAHEDRON  ";
        break;
    case VTK_QUADRATIC_LINEAR_QUAD:
        type = "VTK_QUADRATIC_LINEAR_QUAD ";
        break;
    case VTK_QUADRATIC_LINEAR_WEDGE:
        type = "VTK_QUADRATIC_LINEAR_WEDGE";
        break;
    case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
        type = "VTK_BIQUADRATIC_QUADRATIC_WEDGE ";
        break;
    case VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON:
        type = "VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON";
        break;
    case VTK_BIQUADRATIC_TRIANGLE:
        type = "VTK_BIQUADRATIC_TRIANGLE";
        break;
    case VTK_CUBIC_LINE:
        type = "VTK_CUBIC_LINE";
        break;
    case VTK_CONVEX_POINT_SET:
        type = "VTK_CONVEX_POINT_SET";
        break;
    case VTK_POLYHEDRON:
        type = "VTK_POLYHEDRON";
        break;
    case VTK_PARAMETRIC_CURVE:
        type = "VTK_PARAMETRIC_CURVE";
        break;
    case VTK_PARAMETRIC_SURFACE:
        type = "VTK_PARAMETRIC_SURFACE";
        break;
    case VTK_PARAMETRIC_TRI_SURFACE:
        type = "VTK_PARAMETRIC_TRI_SURFACE";
        break;
    case VTK_PARAMETRIC_QUAD_SURFACE:
        type = "VTK_PARAMETRIC_QUAD_SURFACE";
        break;
    case VTK_PARAMETRIC_TETRA_REGION:
        type = "VTK_PARAMETRIC_TETRA_REGION";
        break;
    case VTK_PARAMETRIC_HEX_REGION:
        type = "VTK_PARAMETRIC_HEX_REGION";
        break;
    case VTK_HIGHER_ORDER_EDGE:
        type = "VTK_HIGHER_ORDER_EDGE";
        break;
    case VTK_HIGHER_ORDER_TRIANGLE:
        type = "VTK_HIGHER_ORDER_TRIANGLE";
        break;
    case VTK_HIGHER_ORDER_QUAD:
        type = "VTK_HIGHER_ORDER_QUAD";
        break;
    case VTK_HIGHER_ORDER_POLYGON:
        type = "VTK_HIGHER_ORDER_POLYGON";
        break;
    case VTK_HIGHER_ORDER_TETRAHEDRON:
        type = "VTK_HIGHER_ORDER_TETRAHEDRON";
        break;
    case VTK_HIGHER_ORDER_WEDGE:
        type = "VTK_HIGHER_ORDER_WEDGE";
        break;
    case VTK_HIGHER_ORDER_PYRAMID:
        type = "VTK_HIGHER_ORDER_PYRAMID";
        break;
    case VTK_HIGHER_ORDER_HEXAHEDRON:
        type = "VTK_HIGHER_ORDER_HEXAHEDRON";
        break;
    case VTK_LAGRANGE_CURVE:
        type = "VTK_LAGRANGE_CURVE";
        break;
    case VTK_LAGRANGE_TRIANGLE:
        type = "VTK_LAGRANGE_TRIANGLE ";
        break;
    case VTK_LAGRANGE_QUADRILATERAL:
        type = "VTK_LAGRANGE_QUADRILATERAL";
        break;
    case VTK_LAGRANGE_TETRAHEDRON:
        type = "VTK_LAGRANGE_TETRAHEDRON";
        break;
    case VTK_LAGRANGE_HEXAHEDRON:
        type = "VTK_LAGRANGE_HEXAHEDRON";
        break;
    case VTK_LAGRANGE_WEDGE:
        type = "VTK_LAGRANGE_WEDGE";
        break;
    case VTK_LAGRANGE_PYRAMID:
        type = "VTK_LAGRANGE_PYRAMID";
        break;
    default:
        type = "UNKNOWN";
        break;
    }

    return type;
}
