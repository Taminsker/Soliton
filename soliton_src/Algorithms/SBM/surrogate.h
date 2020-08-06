#ifndef SURROGATE_H
#define SURROGATE_H

#include <Core/Defs4Soliton/defs4soliton.h>

typedef enum
{
    OUTSIDE,
    INSIDE,
    MIXED
} LEVEL_SET_TAG;

class Mesh;
class Sto4Sol;

SOLITON_RETURN AddLevelSetAndTag (Mesh* mesh, Mesh* object);

SOLITON_RETURN BuildSurrogateDomains (Sto4Sol* store);

#endif // SURROGATE_H
