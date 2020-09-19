#ifndef SBM_H
#define SBM_H

#include <ProgDef/proddef.h>
#include <IO/parsers.h>

enum class TAG_INTERSECTION
{
    TAG_UNKNOWN     = 0,
    TAG_INSIDE      = 1,
    TAG_OUTSIDE     = 2,
    TAG_MIXED       = 3,
    FIRST           = TAG_UNKNOWN,
    LAST            = TAG_MIXED
};


std::string ToString (TAG_INTERSECTION tag);
std::ostream& operator<< (std::ostream& out, TAG_INTERSECTION tag);

class Mesh;
class Point;

void AddLevelSetBetween (Mesh* mesh, Mesh* object);

void TagCellsFromLevelSet (Mesh* mesh, Mesh* object);

void TagEdgesFromTagCells (Mesh* mesh, Mesh* object);


[[deprecated("Use BuildDisplacementVectorsBounds() instead.")]]
void BuildDisplacementVectorsBounds2 (Mesh* mesh, Mesh* object);
void BuildDisplacementVectorsBounds (Mesh* mesh, Mesh* object);

int GetDisplacementVectorAtPoint (Point* atpoint, Mesh* targetToMap, Point* d);

#endif // SBM_H
