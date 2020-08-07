#ifndef HASH4EDGES_H
#define HASH4EDGES_H

#include <algorithm>
#include <unordered_map>
#include <vector>

#include <Core/Defs4Soliton/defs4soliton.h>

class Mesh;
class Cell;
class Point;

typedef struct {
    std::string id = "";
    Cell* cell = nullptr;
    std::vector <Point*> listpoints = {};
    unsigned int gmshcelltype = 0;
} SolitonHashCell;

typedef std::unordered_multimap <std::string, SolitonHashCell*> SolitonHashMap;

std::vector <SolitonHashCell *> ExtractUndergroundCells (Cell* cell);

SOLITON_RETURN AddToHashMap (SolitonHashMap* map, std::vector <SolitonHashCell *>* HashCells);

SOLITON_RETURN BuildEdgesWithHashMap (Mesh* mesh);

//long int HashFunction (std::vector <int> idx, int KeyGenerator);
std::string HashFunction (std::vector <Point*>* pointlist);

void Print (SolitonHashMap* map, std::string name);

#endif // HASH4EDGES_H
