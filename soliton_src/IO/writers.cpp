#include <fstream>

#include <Core/core.h>
#include "writers.h"

SOLITON_RETURN WriterVTK::WithCells (Mesh* mesh)
{
    HEADERFUN("WriterVTK::WithCells");

    DataContainer<Point*>* pointsData = mesh->GetPointsData ();
    DataContainer<Cell*>* cellsData = mesh->GetCellsData ();

    std::string filename = mesh->GetName ()+"_wcells.vtk";

    BEGIN << "Write the vtk file with cells p.o.v : " << filename << ENDLINE;

    std::ofstream outfile (filename);
    int numpoints = mesh->GetNumberOfPoints ();
    int numcells = mesh->GetNumberOfCells ();
    int numInfosCells = mesh->GetNumberOfInfosCells ();

    if (!outfile.is_open ())
    {
        ERROR << "the file \"" << filename << "\" can not be open." << BLINKRETURN << ENDLINE;
        return SOLITON_FAILURE;
    }
    else
    {
#ifdef VERBOSE
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
    INFOS << "written : " << COLOR_YELLOW;
    INFOS << "points ";
#endif

    outfile << "CELLS " << numcells << " " << numInfosCells << std::endl;

    // CELLS
    for (int i = 0; i < numcells; ++i)
        outfile << *mesh->GetCell (i);
    outfile << std::endl;

#ifdef VERBOSE
    INFOS << "cells ";
#endif

    outfile << "CELL_TYPES " << numcells << std::endl;

    // CELLSTYPE
    for (int i = 0; i < numcells; ++i)
        outfile << mesh->GetCell (i)->GetTypeVTK () << std::endl;
    outfile << std::endl;

#ifdef VERBOSE
    INFOS << "cells type " << ENDLINE;
#endif

    // POINTDATA
    outfile << "POINT_DATA " << numpoints << std::endl;

#ifdef VERBOSE
    INFOS << "point data : " << COLOR_YELLOW;
#endif

    for (auto arr : *pointsData->GetIntArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *pointsData->GetBooleanArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " bool " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *pointsData->GetDoubleArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *pointsData->GetVecArrays ()->GetAll ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }
#ifdef VERBOSE
    std::cout << ENDLINE;
#endif

    // CELLDATA
    outfile << "CELL_DATA " << numcells << std::endl;

#ifdef VERBOSE
    INFOS << "cell data : " << COLOR_YELLOW;
#endif
    for (auto arr : *cellsData->GetIntArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *cellsData->GetBooleanArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " bool " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *cellsData->GetDoubleArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *cellsData->GetVecArrays ()->GetAll ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    outfile.close ();

#ifdef VERBOSE
    std::cout << ENDLINE;
    ENDFUN;
#endif

    return SOLITON_SUCCESS;
}

SOLITON_RETURN WriterVTK::WithEdges (Mesh* mesh)
{
    HEADERFUN("WriterVTK::WithEdges");

    DataContainer<Point *>* pointsData = mesh->GetPointsData ();
    DataContainer<Edge *>* edgesData = mesh->GetEdgesData ();

    std::string filename = mesh->GetName ()+"_wedges.vtk";

    BEGIN << "Write the vtk file with edges p.o.v : " << filename << ENDLINE;

    std::ofstream outfile (filename);
    int numpoints = mesh->GetNumberOfPoints ();
    int numedges = mesh->GetNumberOfEdges ();
    int numInfosEdges = mesh->GetNumberOfInfosEdges ();

    if (!outfile.is_open ())
    {
        ERROR << "the file \"" << filename << "\" can not be open." << BLINKRETURN << ENDLINE;
        return SOLITON_FAILURE;
    }
    else
    {
#ifdef VERBOSE
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
    INFOS << "written : " << COLOR_YELLOW;
    INFOS << "points ";
#endif

    outfile << "CELLS " << numedges << " " << numInfosEdges << std::endl;

    // CELLS
    for (int i = 0; i < numedges; ++i)
        outfile << *mesh->GetEdge (i);
    outfile << std::endl;

#ifdef VERBOSE
    INFOS << "edges ";
#endif

    outfile << "CELL_TYPES " << numedges << std::endl;

    // CELLSTYPE
    for (int i = 0; i < numedges; ++i)
        outfile << mesh->GetEdge (i)->GetTypeVTK () << std::endl;
    outfile << std::endl;

#ifdef VERBOSE
    INFOS << "edges type " << ENDLINE;
#endif

    // POINTDATA
    outfile << "POINT_DATA " << numpoints << std::endl;

#ifdef VERBOSE
    INFOS << "point data : " << COLOR_YELLOW;
#endif

    for (auto arr : *pointsData->GetIntArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *pointsData->GetBooleanArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *pointsData->GetDoubleArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *pointsData->GetVecArrays ()->GetAll ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }
    std::cout << ENDLINE;

    // CELLDATA
    outfile << "CELL_DATA " << numedges << std::endl;

#ifdef VERBOSE
    INFOS << "edge data : " << COLOR_YELLOW;
#endif

    for (auto arr : *edgesData->GetIntArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *edgesData->GetBooleanArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " int " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *edgesData->GetDoubleArrays ()->GetAll ())
    {
        outfile << "SCALARS " << arr->name << " double " << std::endl;
        outfile << "LOOKUP_TABLE default" << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

    for (auto arr : *edgesData->GetVecArrays ()->GetAll ())
    {
        outfile << "VECTORS " << arr->name << " double " << std::endl;
        outfile << arr->vec << std::endl;
        outfile << std::endl;
#ifdef VERBOSE
        INFOS << arr->name << " ";
#endif
    }

#ifdef VERBOSE
    std::cout << ENDLINE;
    ENDFUN;
#endif

    outfile.close ();

    return SOLITON_SUCCESS;
}
