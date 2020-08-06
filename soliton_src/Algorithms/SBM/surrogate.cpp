#include "surrogate.h"

#include "../AlgoMesh/algomesh.h"

#include <Core/core.h>
#include <IO/io.h>
#include <Solver/solver.h>

SOLITON_RETURN BuildSurrogateDomains (Sto4Sol* store)
{
    HEADERFUN("BuildSurrogateDomains");
    BEGIN << "Build surrogate domains from all the objects" << ENDLINE;

    Mesh* mesh = store->mesh;
    std::string namemesh = mesh->GetName ();
    int numPointsMesh = mesh->GetNumberOfPoints ();
    //    int numCellsMesh = mesh->GetNumberOfCells ();
    int numEdgesMesh = mesh->GetNumberOfEdges ();
    int objCount = int(store->listobjects.size ());

    // Loop on objects (other meshes)
    for (int objId = 0; objId < objCount; ++objId)
    {
        Mesh* object = store->listobjects.at (std::size_t (objId));
        std::string nameobject = object->GetName ();

        int numCellsObject = object->GetNumberOfCells ();

        auto tagEdgeSurrogateVec = mesh->GetEdgesData ()->GetIntArrays ()->Get (nameobject + NAME_SURROGATETAG);

        if (tagEdgeSurrogateVec == nullptr)
        {
            ERROR << "mesh : you need to tag the edge before use tags to build a surrogate domain..." << BLINKRETURN << ENDLINE;
            return SOLITON_FAILURE;
        }

        auto tagCellVec = object->GetCellsData ()->GetIntArrays ()->Get (namemesh + NAME_INTERSECTIONTAG);

        if (tagCellVec == nullptr)
        {
            ERROR << "object : you need to tag the cell before use tags to build a surrogate domain..." << BLINKRETURN << ENDLINE;
            return SOLITON_FAILURE;
        }

        auto normalCell = object->GetCellsData ()->GetVecArrays ()->Get (NAME_NORMALONCELLS);

        if (normalCell == nullptr)
        {
            ERROR << "object : you need to compute normals on cells before use it..." << BLINKRETURN << ENDLINE;
            return SOLITON_FAILURE;
        }

        // ************************************************************* //

        // Initialize new vec

        std::vector <Point*> d_vec;
        d_vec.resize (std::size_t (numPointsMesh));
        for (int i = 0; i < numPointsMesh; ++i)
            d_vec.at (std::size_t (i)) = new Point();

        std::vector <bool> alreadyTreated (std::size_t (numPointsMesh), false);

        // ************************************************************* //

        int count_error = 0;
        int count_d_vec = 0;

        // Loop on edges of the mesh
        for (int edgeId = 0; edgeId < numEdgesMesh; ++edgeId)
        {
            Edge* edge = mesh->GetEdge (edgeId);

            // Test if the edges is on mixed cell
            if (tagEdgeSurrogateVec->vec.at (std::size_t (edgeId)) == MIXED)
            {
                // Compute d vectors
                // if we are now, the current edge is a surrogate edge ! We need
                // to find the 'best' cell on physical domain to map all points
                // of the edge.

                std::vector <Point*>* listPointsOnEdge = edge->GetPoints ();

                for (Point* point : *listPointsOnEdge)
                {
                    int pointId = point->GetGlobalIndex ();

                    // this point is already treated, so skip it
                    if (alreadyTreated.at (std::size_t (pointId)))
                        continue;

                    std::vector <Point*>* listPointsOnBestCell = nullptr;
                    Point* normal = nullptr;

                    double distMinSum = 1e6;

                    for (int cellId = 0; cellId < numCellsObject; ++cellId)
                    {
                        if (tagCellVec->vec.at (std::size_t (cellId)) == MIXED)
                        {
                            Cell* cell = object->GetCell (cellId);

                            double distTemp = 0;

                            std::vector <Point*>* listPointsOnCell = cell->GetPoints ();

                            for (Point* pointOfCell : *listPointsOnCell)
                                distTemp += EuclidianDist (*pointOfCell, *point);

                            if (distTemp < distMinSum)
                            {
                                distMinSum = distTemp;
                                listPointsOnBestCell = listPointsOnCell;
                                normal = normalCell->vec.at (std::size_t (cellId));
                            }
                        }
                    }

                    // d =  x - \tilde{x}
                    // \tilde{x} = P_0 + sum_{i=1}^{N_b} <P_0 - x|P_0 - P_i> / <P_0 - P_i|P_0 - P_i> (P_0 - P_i)

                    Point* P_0 = listPointsOnBestCell->at (0);
                    Point tildex = *P_0;
                    Point AB = *point - *P_0;

                    for (Point* P_i : *listPointsOnBestCell)
                    {
                        if (P_0 != P_i)
                        {
                            Point V = *P_0 - *P_i;
                            tildex = tildex +  (AB | V) / (V | V) * V;
                        }
                    }

                    Point* d = d_vec.at (std::size_t (point->GetGlobalIndex ()));
                    *d = tildex - *point;
                    *d = (*d | *normal) * *normal;

                    if (std::isnan(d->x) || std::isnan(d->z) || std::isnan(d->z))
                        count_error++;

                    count_d_vec++;
                    alreadyTreated.at (std::size_t (point->GetGlobalIndex ())) = true;
                }
            }
        }
        if (count_error != 0)
            ERROR << "there is " << count_error << " \'nan\' in d vectors..." << ENDLINE;

        INFOS << "displacement vector d for points on surrogate boundary : " << count_d_vec  << ENDLINE;

        mesh->GetPointsData ()->GetVecArrays ()->Add (nameobject + NAME_DISPLACEMENTVECTOR, d_vec);
    }

    ENDFUN;
    return SOLITON_SUCCESS;
}


SOLITON_RETURN AddLevelSetAndTag (Mesh* mesh, Mesh* object)
{
    HEADERFUN("AddLevelSetAndTag");
    BEGIN << "Add level set distance to object and tag in/out : " << COLOR_BLUE << mesh->GetName () << " <-> " << object->GetName () << ENDLINE;

    int numEdgesMesh = mesh->GetNumberOfEdges ();
    int numCellsMesh = mesh->GetNumberOfCells ();
    int numPointsMesh = mesh->GetNumberOfPoints ();
    int numPointsObject = object->GetNumberOfPoints ();
    std::string nameobject = object->GetName ();
    int count = 0;

    std::vector <double> levelSetPoints (std::size_t (numPointsMesh), 0);

    auto normals = object->GetPointsData ()->GetVecArrays ()->Get ("NormalOnPoints");
    if (normals == nullptr)
    {
        ERROR << "you need ton compute normal on point, before use it... " << BLINKRETURN << ENDLINE;
        return SOLITON_FAILURE;
    }

    // levelSetPoints

    for (int i = 0; i < numPointsMesh; ++i)
    {
        Point* p_cur = mesh->GetPoint (i);

        Point* p_min = object->GetPoint (0);
        Point* normal = normals->vec.at (0);

        double dist = EuclidianDist (*p_min, *p_cur);

        for (int j = 0; j < numPointsObject; ++j)
        {

            Point* p_curv = object->GetPoint (j);

            double d_temp = EuclidianDist (*p_curv, *p_cur);

            if (d_temp < dist)
            {
                normal = normals->vec.at (std::size_t (j));
                p_min = p_curv;
                dist = d_temp;
            }
        }

        levelSetPoints.at (std::size_t (i)) = (*normal | (*p_cur - *p_min));
    }

    mesh->GetPointsData ()->GetDoubleArrays ()->Add (nameobject + NAME_LEVELSET, levelSetPoints);

#ifdef VERBOSE
    INFOS << "level-set computed." << ENDLINE;
#endif

    // tagCellInOut
    std::vector <int> tagCellInOut (std::size_t (numCellsMesh), MIXED);

    for (int i = 0; i < numCellsMesh; ++i)
    {
        Cell* c = mesh->GetCell (i);

        int countpstv = 0;
        int countneg = 0;

        for (Point* p : *c->GetPoints ())
        {
            int idx = p->GetGlobalIndex ();
            double d = levelSetPoints.at (std::size_t (idx));
            if (d <= 0)
                countneg++;
            else
                countpstv++;
        }

        if (countneg == 0 && countpstv > 0)
            tagCellInOut.at (std::size_t(i)) = OUTSIDE;
        else if (countpstv == 0 && countneg > 0)
            tagCellInOut.at (std::size_t(i)) = INSIDE;
        else
        {
            tagCellInOut.at (std::size_t(i)) = MIXED;
            count++;
        }
    }

    if (count != 0)
    {
        INFOS << "cell tagging is done : " << count << " cells on object border." << ENDLINE;
    }
    else
    {
        ERROR << REVERSE << "the tagging of the surrogate cells is done but it seems that there is no intersection: the continuation is uncertain but we mark alls the cells to continue." << ENDLINE;

        for (int id = 0; id < numCellsMesh; ++id)
            tagCellInOut.at (std::size_t (id)) = MIXED;
    }

    mesh->GetCellsData ()->GetIntArrays ()->Add (nameobject + NAME_INTERSECTIONTAG, tagCellInOut);


    // tagEdgeSurrogate

    std::vector <int> tagEdgeSurrogate (std::size_t (numEdgesMesh), MIXED);
    count = 0;

    for (int edgeId = 0; edgeId < numEdgesMesh; ++edgeId)
    {
        Edge* edge = mesh->GetEdge (edgeId);
        int* value = &tagEdgeSurrogate.at (std::size_t (edgeId));

        std::vector <Cell*>* cellList = edge->GetCellsList ();

        int count_INSIDE = 0;
        int count_OUTSIDE = 0;
        int count_MIXED = 0;

        for (Cell* cell : *cellList)
            switch (tagCellInOut.at (std::size_t (cell->GetGlobalIndex ())))
            {
            case INSIDE:
                count_INSIDE++;
                break;
            case MIXED:
                count_MIXED++;
                break;
            default:
                count_OUTSIDE++;
                break;
            }


        if (count_MIXED != 0)
        {
            if (count_INSIDE != 0)
            {

                *value = MIXED;
                count++;
            }
            else
                *value = OUTSIDE;
        }
        else
        {
            if (count_INSIDE != 0)
                *value = INSIDE;
            else
                *value = OUTSIDE;
        }
    }



    if (count != 0)
        INFOS << "surrogate edge tagging is done : " << count << " surrogate edges" << ENDLINE;
    else
    {
        ERROR << REVERSE << "the tagging of the surrogate edges is done but it seems that there is no intersection: the continuation is uncertain but we mark alls the edges to continue." << ENDLINE;

        for (int id = 0; id < numEdgesMesh; ++id)
            tagEdgeSurrogate.at (std::size_t (id)) = MIXED;
    }

    mesh->GetEdgesData ()->GetIntArrays ()->Add (nameobject + NAME_SURROGATETAG, tagEdgeSurrogate);

    ENDFUN;
    return SOLITON_SUCCESS;

}
