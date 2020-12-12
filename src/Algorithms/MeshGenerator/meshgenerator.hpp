#ifndef SRC_ALGORITHMS_MESHGENERATOR_MESHGENERATOR_HPP
#define SRC_ALGORITHMS_MESHGENERATOR_MESHGENERATOR_HPP

#include <string>

#include "../../solitonheader.hpp"

class MeshStorage;
class InputDatStruct;
class Mesh;
struct ObjectDatStruct;

std::string GenerateWithGMSH (InputDatStruct * input);

void ObjectGenerator (InputDatStruct * data, MeshStorage * out);

void AlgoGen1 (ObjectDatStruct * data, Mesh * object);
void AlgoGen2 (ObjectDatStruct * data, Mesh * object);
void AlgoGen3 (ObjectDatStruct * data, Mesh * object);

#endif /* SRC_ALGORITHMS_MESHGENERATOR_MESHGENERATOR_HPP */
