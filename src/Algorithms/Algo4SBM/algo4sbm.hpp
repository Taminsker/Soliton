#ifndef SRC_ALGORITHMS_ALGO4SBM_ALGO4SBM_HPP
#define SRC_ALGORITHMS_ALGO4SBM_ALGO4SBM_HPP

#include "../../solitonheader.hpp"

enum class INTER : ul_t;
class Mesh;
class Point;

void AddLevelSetBetween (Mesh * mesh, Mesh * object);

void TagCellsFromLevelSet (Mesh * mesh, Mesh * object);

void TagEdgesFromTagCells (Mesh * mesh, Mesh * object);

[[deprecated ("Use BuildDisplacementVectorsBounds() instead.")]] void BuildDisplacementVectorsBounds2 (Mesh * mesh, Mesh * object);

void BuildDisplacementVectorsBounds (Mesh * mesh, Mesh * object);

int GetDisplacementVectorAtPoint (Point * atpoint, Mesh * targetToMap, Point * d);

#endif /* SRC_ALGORITHMS_ALGO4SBM_ALGO4SBM_HPP */
