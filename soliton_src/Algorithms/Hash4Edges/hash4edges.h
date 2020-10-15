#ifndef HASH4EDGES_H
#define HASH4EDGES_H

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <ProgDef/proddef.h>
#include <Enums/enums.h>
class Mesh;
class Cell;
class Point;
class Edge;

typedef std::unordered_multimap <std::string, Edge*> SolitonHashMap;


void BuildEdgesWithHashMap (Mesh* mesh);
std::string HashFunction (std::vector<int>* idx);
void Print (SolitonHashMap* map, std::string name);

#endif // HASH4EDGES_H
