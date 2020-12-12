#ifndef SRC_IO_WRITERS_HPP
#define SRC_IO_WRITERS_HPP

#include <string>

#include "../solitonheader.hpp"

class Mesh;

void WriteVTKWithEdges (Mesh * mesh, std::string add2basename = "", bool view = true);
void WriteVTKWithCells (Mesh * mesh, std::string add2basename = "", bool view = true);

#endif /* SRC_IO_WRITERS_HPP */
