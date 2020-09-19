#ifndef HASH4EDGES_H
#define HASH4EDGES_H

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <ProgDef/proddef.h>
#include <Core/Cell/celltype.h>

class Mesh;
class Cell;
class Point;

typedef struct {
    std::string id = "";
    Cell* cell = nullptr;
    std::vector <Point*> listpoints = {};
    GMSH_CELL_TYPE gmshcelltype = GMSH_CELL_TYPE::GMSH_0_NODE_EMPTY;
} SolitonHashCell;

typedef std::unordered_multimap <std::string, SolitonHashCell*> SolitonHashMap;

std::vector <SolitonHashCell *> ExtractUndergroundCells (Cell* cell);

void AddToHashMap (SolitonHashMap* map, std::vector <SolitonHashCell *>* HashCells);

void BuildEdgesWithHashMap (Mesh* mesh);

//long int HashFunction (std::vector <int> idx, int KeyGenerator);
std::string HashFunction (std::vector <Point*>* pointlist);

void Print (SolitonHashMap* map, std::string name);

#endif // HASH4EDGES_H
