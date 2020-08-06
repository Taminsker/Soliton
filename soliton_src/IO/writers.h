#ifndef WRITERS_H
#define WRITERS_H

#include <string>
#include <Core/Defs4Soliton/defs4soliton.h>

class Mesh;

namespace WriterVTK {
    SOLITON_RETURN WithEdges (Mesh* mesh);
    SOLITON_RETURN WithCells (Mesh* mesh);
}


#endif // WRITERS_H
