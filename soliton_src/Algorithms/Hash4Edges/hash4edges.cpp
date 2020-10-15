#include "hash4edges.h"

#include <Core/core.h>
#include <Solver/FEStruct/festruct.h>
#include <fstream>

#include <algorithm>


std::string HashFunction (std::vector<int>* idx)
{
    std::string out = "";

    std::sort (idx->begin (), idx->end ());

    for (std::size_t id = 0; id < idx->size (); ++id)
        out += std::to_string (idx->at(id))+",";

    return out;
}

void BuildEdgesWithHashMap (Mesh* mesh)
{
    HEADERFUN("BuildEdgesWithHashMap");

    BEGIN << "Build edges with a hash map on the mesh : " << COLOR_BLUE << mesh->GetName () << ENDLINE;

    SolitonHashMap map;
    FEStore store;
    int numCells = mesh->GetNumberOfCells ();
    int globalIdx = 0;

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);
        std::vector<Point*>* ptsOnCell = cell->GetPoints ();

        FEBase* febase = store.GetElementFor (cell);

        for (std::size_t edgeId = 0; edgeId < febase->GetNumberOfEdges (); ++edgeId)
        {
            Edge* edge = febase->GetEdge (edgeId);

            std::vector<int>    idxPoints = {};
            std::vector<Point*> listPoints = {};
            for (Point* pt : *edge->GetPoints ())
            {
                listPoints.push_back (ptsOnCell->at (static_cast<std::size_t>(pt->GetGlobalIndex ())));
                idxPoints.push_back (listPoints.back ()->GetGlobalIndex ());
            }

            std::string key = HashFunction (&idxPoints);

            SolitonHashMap::iterator list = map.find (key);

            if (list != map.end ())
                list->second->AddCell (cell);
            else
            {
                Edge* newEdge = new Edge();
                newEdge->AddPoints (listPoints);
                newEdge->SetTag (edge->GetTag ());
                newEdge->SetType (edge->GetTypeGMSH ());
                newEdge->SetGlobalIndex (globalIdx);
                newEdge->AddCell (cell);

                map.insert ({key, newEdge});
                globalIdx++;
            }
        }
    }

#ifdef PRINTHASHMAP
    Print(&map, "hashmap_"+mesh->GetName()+".txt");
#endif

#ifdef VERBOSE
    INFOS << "Hash table size : " << map.size () << ENDLINE;
#endif

    for (std::pair<std::string, Edge*> obj : map)
        mesh->AddEdge (obj.second);

    INFOS << "Build : " << mesh->GetNumberOfEdges () << " edges on mesh " << COLOR_BLUE << "\"" << mesh->GetName () << COLOR_DEFAULT << "\"." << ENDLINE;

    map.clear ();

    ENDFUN;
    return;
}

void Print (SolitonHashMap* map, std::string name)
{
    VOID_USE(map); VOID_USE(name);

    HEADERFUN ("Print HashMap");

    std::ofstream output (name);

    if (!output.is_open ())
        INFOS << "Can not be open" << ENDLINE;

    for (auto obj : *map)
    {
        Edge* edge = obj.second;

        output << "\tID: " << obj.first << "        \t idx-point " << std::flush;

        for (Point* point : *edge->GetPoints ())
            output << point->GetGlobalIndex () << ", " << std::flush;

        output << "         \tGLOBALIDX-Cell : " << std::flush;
        for (Cell* cell : *edge->GetCellsList ())
            output << cell->GetGlobalIndex () << ", " << std::flush;

        output << "         \tGMSHTYPE-Cell : " << std::flush;
        for (Cell* cell : *edge->GetCellsList ())
            output << to_string(cell->GetTypeGMSH ()) << ", " << std::flush;

        output << "         \tGMSHTYPE-EDGE : " << std::flush;
        output << to_string(edge->GetTypeGMSH ()) << '\t' << std::endl;
        output << "NEXT HASH VECTOR "<< SEPARATOR << std::endl;

    }

    INFOS << "Print hash map in \"" << COLOR_BLUE << name << COLOR_DEFAULT << "\"" << ENDLINE;

    output.close ();

    return;
}
