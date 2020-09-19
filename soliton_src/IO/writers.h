#ifndef WRITERS_H
#define WRITERS_H

#include <string>
#include <ProgDef/proddef.h>

class Mesh;


void WriteVTKWithEdges (Mesh* mesh, std::string add2basename = "");
void WriteVTKWithCells (Mesh* mesh, std::string add2basename = "");


#endif // WRITERS_H
