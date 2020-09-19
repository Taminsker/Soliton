#ifndef MESHGENERATOR_H
#define MESHGENERATOR_H

#include <string>

#include <ProgDef/proddef.h>

class Sto4Sol;
class InputDatStruct;
class Mesh;
struct ObjectDatStruct;

std::string GenerateWithGMSH(InputDatStruct* input);

void ObjectGenerator (InputDatStruct* data, Sto4Sol* out);

void AlgoGen1 (ObjectDatStruct* data, Mesh* object);
void AlgoGen2 (ObjectDatStruct* data, Mesh* object);
void AlgoGen3 (ObjectDatStruct* data, Mesh* object);

#endif // MESHGENERATOR_H