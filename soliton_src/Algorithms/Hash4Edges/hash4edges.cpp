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
    BEGIN << "Build edges with a hash map on the mesh : " << COLOR_BLUE << mesh->GetName () << ENDLINE;

    SolitonHashMap map;
    FEStore store;
    int numCells = mesh->GetNumberOfCells ();
    int globalIdx = 0;

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);
        FEBase* febase = store.GetElementFor (cell);

        //        ERROR << "Cell id : " << cellId << ENDLINE;

        for (std::size_t edgeId = 0; edgeId < febase->GetNumberOfEdges (); ++edgeId)
        {
            Edge* edge = febase->GetEdge (edgeId);

            //            ERROR << "  edge id-loc : " << edgeId << " with " << edge->GetNumberOfPoints () << " pts" << ENDLINE;


            std::vector<int>    idxPoints = {};
            std::vector<Point*> listPoints = {};

            for (int ptLocId = 0; ptLocId < edge->GetNumberOfPoints (); ++ptLocId)
            {
                Point* ptOnEdge = edge->GetPoints ()->at (static_cast<std::size_t>(ptLocId));
                //                Point* ptOnFebase = febase->GetPoint (static_cast<std::size_t>(ptLocId));
                Point* ptOnCell = cell->GetPoints ()->at (static_cast<std::size_t>(ptOnEdge->GetGlobalIndex ()));

                //                ERROR << "      febase : ptOnEdge : " << ptOnEdge->GetGlobalIndex () << ENDLINE;
                //                ERROR << "      febase : ptGlob : " << ptOnFebase->GetGlobalIndex () << ENDLINE;
                //                ERROR << "      febase : ptOnCell : " << ptOnCell->GetGlobalIndex () << ENDLINE;

                listPoints.push_back (ptOnCell);
                idxPoints.push_back (ptOnCell->GetGlobalIndex ());
            }

            std::string key = HashFunction (&idxPoints);

            //            ERROR << "  key : " << key << ENDLINE;
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

            //            ERROR << ENDLINE;
        }
    }

#ifdef PRINTHASHMAP
    Print(&map, "hashmap_"+mesh->GetName()+".txt");
#endif

#ifdef VERBOSE
    INFOS << "Hash table size : " << map.size () << ENDLINE;
#endif

    int globalIdEdge = 0;
    for (std::pair<std::string, Edge*> obj : map)
    {
        Edge* edge = obj.second;
        edge->SetGlobalIndex (globalIdEdge);
        mesh->AddEdge (edge);

        globalIdEdge++;
    }

    INFOS << "Build : " << mesh->GetNumberOfEdges () << " edges on mesh " << COLOR_BLUE << "\"" << mesh->GetName () << COLOR_DEFAULT << "\"." << ENDLINE;

    map.clear ();

    //    exit(0);
    ENDFUN;
    return;
}

void Print (SolitonHashMap* map, std::string name)
{
    VOID_USE(map); VOID_USE(name);

    std::ofstream output (name);

    if (!output.is_open ())
        INFOS << "Can not be open" << ENDLINE;

    for (auto obj : *map)
    {
        Edge* edge = obj.second;

        for (Cell* cell : *edge->GetCellsList ())
        {
            output << "\tkey: " << obj.first << "edge id : " << edge->GetGlobalIndex () << "        \t idx-pts " << std::flush;

            for (Point* point : *edge->GetPoints ())
                output << point->GetGlobalIndex () << ", " << std::flush;

            output << "         \tGMSHTYPE-EDGE : " << std::flush;
            output << to_string(edge->GetTypeGMSH ()) << '\t' << std::flush;

            output << "         \tGLOBALIDX-Cell : " << std::flush;
            output << cell->GetGlobalIndex () << " " << std::flush;
            output << "         \tGMSHTYPE-Cell : " << std::flush;
            output << to_string(cell->GetTypeGMSH ()) << " " << std::endl;

        }
        output << "NEXT HASH VECTOR "<< SEPARATOR << std::endl;

    }

    INFOS << "Print hash map in \"" << COLOR_BLUE << name << COLOR_DEFAULT << "\"" << ENDLINE;

    output.close ();

    return;
}
