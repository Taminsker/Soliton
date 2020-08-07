#include "hash4edges.h"

#include <Core/core.h>
#include <fstream>

#include <algorithm>


//long int HashFunction (std::vector <Point*> pointlist, int KeyGenerator)
//{
//    std::size_t numidx = pointlist.size ();
//    long int out = 0;

//    std::sort (pointlist.begin (), pointlist.end ());

//    for (std::size_t id = 0; id < numidx; ++id)
//        out += static_cast<long int>(KeyGenerator)^static_cast<long int>(id + 1) * static_cast<long int>(pointlist.at (id)->GetGlobalIndex ());


//    return out;

//}

std::string HashFunction (std::vector <Point*>* pointlist)
{
    std::size_t numidx = pointlist->size ();
    std::string out = "";

    std::sort (pointlist->begin (), pointlist->end ());

    for (std::size_t id = 0; id < numidx; ++id)
        out += std::to_string (pointlist->at (id)->GetGlobalIndex ())+",";


    return out;

}

std::vector <SolitonHashCell *> ExtractUndergroundCells (Cell* cell)
{
    HEADERFUN("ExtractUndergroundCells");

    std::vector <SolitonHashCell *> out;
    std::size_t np = std::size_t (cell->GetNumberOfPoints ());
    auto lp = cell->GetPoints ();
    unsigned int type = cell->GetTypeVTK ();

    switch (type)
    {
    case VTK_LINE:
    {
        unsigned int edgeType = 15; // Vertex

        for (std::size_t j = 0; j < np; j++)
        {
            SolitonHashCell* add = new SolitonHashCell();

            add->gmshcelltype = edgeType;
            add->cell = cell;
            add->listpoints.push_back (lp->at (j));
            out.push_back (add);
        }
        break;
    }
    case VTK_QUADRATIC_EDGE:
    {
        unsigned int edgeType = 15; // Vertex

        for (std::size_t j = 0; j < np-2; j=j+2)
        {
            SolitonHashCell* add = new SolitonHashCell();

            add->gmshcelltype = edgeType;
            add->cell = cell;
            add->listpoints.push_back (lp->at (j));
            out.push_back (add);
        }
        break;
    }
    case VTK_TRIANGLE:
    {
        unsigned int edgeType = 1; // LINE

        for (std::size_t j = 0; j < np; ++j)
        {
            SolitonHashCell* add = new SolitonHashCell();

            add->gmshcelltype = edgeType;
            add->cell = cell;

            if (j != np -1)
                add->listpoints = {lp->at (j), lp->at (j+1)};
            else
                add->listpoints = {lp->at (j), lp->at (0)};

            out.push_back (add);
        }

        break;
    }
    case VTK_QUADRATIC_TRIANGLE:
    {
        unsigned int edgeType = 8; // second order line

        for (std::size_t j = 0; j < np - 1; ++j)
        {
            SolitonHashCell* add = new SolitonHashCell();

            add->gmshcelltype = edgeType;
            add->cell = cell;

            if (j != np -2)
                add->listpoints = {lp->at (j), lp->at (j+1), lp->at (j+2)};
            else
                add->listpoints = {lp->at (j), lp->at (j+1), lp->at (0)};

            out.push_back (add);
        }

        break;
    }
    default:
    {
        ERROR << "the cell (type " << type << ") is not supported yet.. please look the meshtools file in ExtractBoundary." << ENDLINE;
        break;
    }
    }

    return out;
}

SOLITON_RETURN AddToHashMap (SolitonHashMap* map, std::vector <SolitonHashCell *>* HashCells)
{
    HEADERFUN("AddToHashMap");

    for (SolitonHashCell* obj : *HashCells)
    {
        std::string id = obj->id;
        bool noid = (obj->id == std::to_string (NONEID));
        bool nocell = (obj->cell == nullptr);
        bool nopoints = (obj->listpoints.size () == 0);

        if (noid || nocell || nopoints)
        {
            ERROR << "seems to be an error when we try add to hash map, but we continue..." << ENDLINE;
            continue;
        }

        map->insert ({id, obj});
    }

    return SOLITON_SUCCESS;
}

SOLITON_RETURN BuildEdgesWithHashMap (Mesh* mesh)
{
    HEADERFUN("BuildEdgesWithHashMap");

    BEGIN << "Build edges with a hash map on the mesh : " << COLOR_BLUE << mesh->GetName () << ENDLINE;

    SolitonHashMap map;
//    int numPoints = mesh->GetNumberOfPoints ();
    int numCells = mesh->GetNumberOfCells ();

//    int keyGenerator = numPoints;

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);
        std::vector <SolitonHashCell *> vectorHash = ExtractUndergroundCells (cell);
        // Compute key hash

        for (SolitonHashCell* obj : vectorHash)
                    obj->id = HashFunction (&obj->listpoints);

        SOLITON_RETURN error = AddToHashMap (&map, &vectorHash);
        USE_SOLITON_RETURN(error);
    }

    Print(&map, "test.txt");

#ifdef VERBOSE
    INFOS << "Hash table size : " << map.size () << ENDLINE;
#endif

    Edge* currentedge = nullptr;
    std::string idx = std::to_string (NONEID);
    int globalId = 0;

    std::vector <int> counter;
    counter.resize (map.bucket_count ());

    for (auto obj : map)
    {
        Cell* cell = obj.second->cell;

        if (idx != obj.first)
        {
            if (idx != std::to_string (NONEID) && currentedge != nullptr)
                mesh->AddEdge (currentedge);

            std::string objIdx = obj.first;
            std::vector<Point*> vecPoints = obj.second->listpoints;

            idx = objIdx;
            currentedge = new Edge ();

            for (Point* point : vecPoints)
                currentedge->AddPoint (point);

            currentedge->SetType (obj.second->gmshcelltype);
            currentedge->SetGlobalIndex (globalId);
            globalId++;
        }

        currentedge->AddCell (cell);

        delete obj.second;
        obj.second = nullptr;
    }


    if (idx != std::to_string (NONEID) && currentedge != nullptr)
        mesh->AddEdge (currentedge);

    INFOS << "build : " << mesh->GetNumberOfEdges () << " edges." << ENDLINE;

    map.clear ();

    ENDFUN;
    return SOLITON_SUCCESS;
}

void Print (SolitonHashMap* map, std::string name)
{

    HEADERFUN ("Print HashMap");
    BEGIN << "display hash map" << ENDLINE;

    std::ofstream output (name);

    if (!output.is_open ())
        INFOS << "Can not be open" << ENDLINE;

    for (auto obj : *map)
    {
        output << "ID: " << obj.first << "         \tID': " << obj.second->id << "         \t idx-point " << std::flush;
        for (Point* point : obj.second->listpoints)
            output << point->GetGlobalIndex () << ", " << std::flush;
        output << "         \tGLOBALIDX-Cell : " << obj.second->cell->GetGlobalIndex () << "         \tGMSHTYPE-Cell : " << std::flush;
        output << obj.second->gmshcelltype << '\t' << std::endl;

    }

    INFOS << "PRINT DONE in " << name << ENDLINE;

    output.close ();

    return;
}
