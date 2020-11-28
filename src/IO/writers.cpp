#include <fstream>

#include <Core/core.h>
#include "writers.h"

void WriteVTKWithCells (Mesh* mesh, std::string add2basename, bool view)
{
    DataContainer<Point*>* pointsData = mesh->GetPointsData ();
    DataContainer<Cell*>* cellsData = mesh->GetCellsData ();
    std::string filename = "";

    if (add2basename == "")
        filename = mesh->GetName ()+"_w_cells.vtk";
    else
        filename = mesh->GetName ()+"_w_cells_" + add2basename + ".vtk";

    if (view)
        BEGIN << "Write the vtk file with cells p.o.v : " << COLOR_BLUE << filename << ENDLINE;

    std::ofstream outfile (filename);
    int numpoints = mesh->GetNumberOfPoints ();
    int numcells = mesh->GetNumberOfCells ();
    int numInfosCells = mesh->GetNumberOfInfosCells ();

    if (!outfile.is_open ())
    {
        ERROR << "the file \"" << filename << "\" can not be open." << BLINKRETURN << ENDLINE;
        return;
    }
    else
    {
#ifdef VERBOSE
        if (view)
            STATUS << "the file \"" << COLOR_BLUE << filename << COLOR_DEFAULT << "\" is open." << ENDLINE;
#endif
    }

    outfile << "# vtk DataFile Version 2.0" << std::endl;
    outfile << filename << ", Soliton output." << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outfile << "POINTS " << numpoints << " double" << std::endl;

    // POINTS
    for (int i = 0; i < numpoints; ++i)
        outfile << *mesh->GetPoint (i) << std::endl;

    outfile << std::endl;

#ifdef VERBOSE
    if (view)
        TREE_BRANCH << "physical : " << COLOR_YELLOW << "points" << ENDLINE;
#endif

    outfile << "CELLS " << numcells << " " << numInfosCells << std::endl;

    // CELLS
    for (int i = 0; i < numcells; ++i)
        outfile << *mesh->GetCell (i);
    outfile << std::endl;

#ifdef VERBOSE
    if (view)
        TREE_BRANCH << "physical : " << COLOR_YELLOW << "cells" << ENDLINE;
#endif

    outfile << "CELL_TYPES " << numcells << std::endl;

    // CELLSTYPE
    for (int i = 0; i < numcells; ++i)
        outfile << static_cast<int>(mesh->GetCell (i)->GetTypeVTK ()) << std::endl;
    outfile << std::endl;

#ifdef VERBOSE
    if (view)
        TREE_BRANCH << "physical : " << COLOR_YELLOW << "cells type" << ENDLINE;
#endif

    // POINTDATA
    outfile << "POINT_DATA " << numpoints << std::endl;

    for (auto arr : *pointsData->GetAll<int> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;

#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << " int " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *pointsData->GetAll<bool> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << " bool " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *pointsData->GetAll<real_t> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << "real_t" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *pointsData->GetAll<Point*> ()->All ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << "Point*" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    // CELLDATA
    outfile << "CELL_DATA " << numcells << std::endl;

    for (auto arr : *cellsData->GetAll<int> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "cell data [" << COLOR_RED << " int " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *cellsData->GetAll<bool> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "cell data [" << COLOR_RED << " bool " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *cellsData->GetAll<real_t> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "cell data [" << COLOR_RED << "real_t" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *cellsData->GetAll<Point*> ()->All ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << std::scientific << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "cell data [" << COLOR_RED << "Point*" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    outfile.close ();

#ifdef VERBOSE
    if (view)
        ENDFUN;
#endif

    return;
}

void WriteVTKWithEdges (Mesh* mesh, std::string add2basename, bool view)
{
    DataContainer<Point *>* pointsData = mesh->GetPointsData ();
    DataContainer<Edge *>* edgesData = mesh->GetEdgesData ();

    std::string filename = "";

    if (add2basename == "")
        filename = mesh->GetName ()+"_w_edges.vtk";
    else
        filename = mesh->GetName ()+"_w_edges_" + add2basename + ".vtk";

    if (view)
        BEGIN << "Write the vtk file with edges p.o.v : " << COLOR_BLUE << filename << ENDLINE;

    std::ofstream outfile (filename);
    int numpoints = mesh->GetNumberOfPoints ();
    int numedges = mesh->GetNumberOfEdges ();
    int numInfosEdges = mesh->GetNumberOfInfosEdges ();

    if (!outfile.is_open ())
    {
        ERROR << "the file \"" << filename << "\" can not be open." << BLINKRETURN << ENDLINE;
        return;
    }
    else
    {
#ifdef VERBOSE
        if (view)
            STATUS << "the file \"" << COLOR_BLUE << filename << COLOR_DEFAULT << "\" is open." << ENDLINE;
#endif
    }

    outfile << "# vtk DataFile Version 2.0" << std::endl;
    outfile << filename << ", Soliton output." << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outfile << "POINTS " << numpoints << " double" << std::endl;

    // POINTS
    for (int i = 0; i < numpoints; ++i)
        outfile << *mesh->GetPoint (i) << std::endl;

    outfile << std::endl;

#ifdef VERBOSE
    if (view)
        TREE_BRANCH << "physical : " << COLOR_YELLOW << "points" << ENDLINE;
#endif

    outfile << "CELLS " << numedges << " " << numInfosEdges << std::endl;

    // CELLS
    for (int i = 0; i < numedges; ++i)
        outfile << *mesh->GetEdge (i);
    outfile << std::endl;

#ifdef VERBOSE
    if (view)
        TREE_BRANCH << "physical : " << COLOR_YELLOW << "edges" << ENDLINE;
#endif

    outfile << "CELL_TYPES " << numedges << std::endl;

    // CELLSTYPE
    for (int i = 0; i < numedges; ++i)
        outfile << static_cast<int>(mesh->GetEdge (i)->GetTypeVTK ()) << std::endl;
    outfile << std::endl;

#ifdef VERBOSE
    if (view)
        TREE_BRANCH << "physical : " << COLOR_YELLOW << "edges type" << ENDLINE;
#endif

    // POINTDATA
    outfile << "POINT_DATA " << numpoints << std::endl;

    for (auto arr : *pointsData->GetAll<int> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << " int " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *pointsData->GetAll<bool> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << " bool " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *pointsData->GetAll<real_t> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << "real_t" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *pointsData->GetAll<Point*> ()->All ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "point data [" << COLOR_RED << "Point*" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    // CELLDATA
    outfile << "CELL_DATA " << numedges << std::endl;

    for (auto arr : *edgesData->GetAll<int> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "edge data [" << COLOR_RED << " int " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *edgesData->GetAll<bool> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "edge data [" << COLOR_RED << " bool " << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *edgesData->GetAll<real_t> ()->All ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "edge data [" << COLOR_RED << "real_t" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

    for (auto arr : *edgesData->GetAll<Point*> ()->All ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        if (view)
            TREE_BRANCH << "edge data [" << COLOR_RED << "Point*" << COLOR_DEFAULT << "]: " << COLOR_YELLOW << arr->name << ENDLINE;
#endif
    }

#ifdef VERBOSE
    if (view)
        ENDFUN;
#endif

    outfile.close ();

    return;
}
