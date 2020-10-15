/** \file algomesh.h */

#ifndef ALGOMESH_H
#define ALGOMESH_H

#include <vector>
#include <string>

#include <ProgDef/proddef.h>
#include <IO/parsers.h>

class Mesh;
class Point;
class Cell;
class Edge;

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Build node-to-node connections.
 * @param mesh the mesh on which to build connections.
 * @return void.
 */
void Build_NtoN (Mesh* mesh);

/**
 * @brief Compute th normals vectors on cells. In particular here, it's the 3D normal to a surface which is built.
 * @param mesh the mesh we are working on
 * @return void
 */
void ComputeNormalsOnCells  (Mesh* mesh);

/**
 * @brief Compute th normals vectors on edges. In particular here, it's the 2D normal to an edge which is built.
 * @param mesh the mesh we are working on
 * @return void
 */
void ComputeNormalsOnEdges  (Mesh* mesh);

/**
 * @brief Compute th normals vectors on points. In particular here, it's the 2D normal to a point which is built (we use the the normals to edges to build them).
 * @param mesh the mesh we are working on
 * @return void
 */
void ComputeNormalsOnPoints (Mesh* mesh);

/**
 * @brief Move a mesh with a geometric transformation.
 * @param mesh the mesh we are working on.
 * @param radius of the transformation.
 * @param center set the new center of the mesh.
 * @return void
 */
void MoveObject (Mesh* mesh, double radius, Point center);


void    ComputeTagPhysical (Mesh* mesh, InputDatStruct* struc);
void    ComputeDampingArea (Mesh* mesh, PHYS tag, double h);
double  GetDampingCoeffFor (Point* atpoint, Mesh* mesh, PHYS tag, double h);

/*! @} End of Doxygen Groups*/

#endif // MESHTOOLS_H
