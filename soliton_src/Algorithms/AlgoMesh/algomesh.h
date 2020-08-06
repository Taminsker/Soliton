#ifndef ALGOMESH_H
#define ALGOMESH_H

#include <vector>
#include <string>

#include <Core/Defs4Soliton/defs4soliton.h>

class Mesh;
class Point;
class Cell;
class Edge;

SOLITON_RETURN Build_NtoN (Mesh* mesh);

SOLITON_RETURN ComputeNormalsOnCells  (Mesh* mesh);
SOLITON_RETURN ComputeNormalsOnEdges  (Mesh* mesh);
SOLITON_RETURN ComputeNormalsOnPoints (Mesh* mesh);

SOLITON_RETURN MoveObject (Mesh* mesh, double radius, Point center);

#endif // MESHTOOLS_H
