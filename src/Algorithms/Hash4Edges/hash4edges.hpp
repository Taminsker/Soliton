#ifndef SRC_ALGORITHMS_HASH4EDGES_HASH4EDGES_HPP
#define SRC_ALGORITHMS_HASH4EDGES_HASH4EDGES_HPP

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "../../solitonheader.hpp"

class Mesh;
class Cell;
class Point;
class Edge;

typedef std::unordered_multimap<std::string, Edge *> SolitonHashMap;

void        BuildEdgesWithHashMap (Mesh * mesh);
std::string HashFunction (std::vector<int> * idx);
void        Print (SolitonHashMap * map, std::string name);

#endif /* SRC_ALGORITHMS_HASH4EDGES_HASH4EDGES_HPP */
