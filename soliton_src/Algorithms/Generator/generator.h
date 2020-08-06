#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>

#include <Core/Defs4Soliton/defs4soliton.h>

class Sto4Sol;
class InputDatStruct;
class Mesh;
struct ObjectDatStruct;

namespace GMSH {

    std::string Generate(InputDatStruct* input);
}

SOLITON_RETURN ObjectGenerator (InputDatStruct* data, Sto4Sol* out);

SOLITON_RETURN AlgoGen1 (ObjectDatStruct* data, Mesh* object);
SOLITON_RETURN AlgoGen2 (ObjectDatStruct* data, Mesh* object);
SOLITON_RETURN AlgoGen3 (ObjectDatStruct* data, Mesh* object);

#endif // GENERATOR_H
