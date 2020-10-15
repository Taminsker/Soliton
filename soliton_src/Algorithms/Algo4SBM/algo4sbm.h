#ifndef SBM_H
#define SBM_H

#include <ProgDef/proddef.h>
#include <IO/parsers.h>
#include <Enums/enums.h>

std::string ToString (INTER tag);
std::ostream& operator<< (std::ostream& out, INTER tag);

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
