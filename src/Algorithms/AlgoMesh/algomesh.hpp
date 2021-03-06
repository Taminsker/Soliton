#ifndef SRC_ALGORITHMS_ALGOMESH_ALGOMESH_HPP
#define SRC_ALGORITHMS_ALGOMESH_ALGOMESH_HPP

/** \file algomesh.hpp */

#include <string>
#include <vector>

#include "../../solitonheader.hpp"

class Mesh;
class Point;
class Cell;
class Edge;
class InputDatStruct;
enum class PHYS : ul_t;

/*!
 * \addtogroup Outils
 * @{
 */

/**
 * @brief Build node-to-node connections.
 * @param mesh the mesh on which to build connections.
 * @return void.
 */
void Build_NtoN (Mesh * mesh);

/**
 * @brief Compute th normals vectors on cells. In particular here, it's the 3D normal to a surface which is built.
 * @param mesh the mesh we are working on
 * @return void
 */
void ComputeNormalsOnCells (Mesh * mesh);

/**
 * @brief Compute th normals vectors on edges. In particular here, it's the 2D normal to an edge which is built.
 * @param mesh the mesh we are working on
 * @return void
 */
void ComputeNormalsOnEdges (Mesh * mesh);

/**
 * @brief Compute th normals vectors on points. In particular here, it's the 2D normal to a point which is built (we use the the normals to edges to build them).
 * @param mesh the mesh we are working on
 * @return void
 */
void ComputeNormalsOnPoints (Mesh * mesh);

/**
 * @brief Move a mesh with a geometric transformation.
 * @param mesh the mesh we are working on.
 * @param radius of the transformation.
 * @param center set the new center of the mesh.
 * @return void
 */
void MoveObject (Mesh * mesh, real_t radius, Point center);

void   ComputeTagPhysical (Mesh * mesh, InputDatStruct * struc);
void   ComputeDampingArea (Mesh * mesh, PHYS tag, real_t h);
real_t GetDampingCoeffFor (Point * atpoint, Mesh * mesh, PHYS tag, real_t h);

/*! @} End of Doxygen Groups*/

#endif /* SRC_ALGORITHMS_ALGOMESH_ALGOMESH_HPP */
