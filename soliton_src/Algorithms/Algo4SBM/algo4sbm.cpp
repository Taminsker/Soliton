#include <Algorithms/Algo4SBM/algo4sbm.h>
#include <Core/core.h>
#include <IO/io.h>
#include <Solver/solver.h>

std::string ToString (INTER tag)
{
    switch (tag)
    {
    case INTER::UNKNOWN:
        return "UNKNOWN";
    case INTER::IN:
        return "IN";
    case INTER::OUT:
        return "OUT";
    case INTER::MIXED:
        return "MIXED";
    }

    return "UNKNOWN";
}

std::ostream& operator<< (std::ostream& out, INTER tag)
{
    out << ToString (tag);
    return out;
}


void BuildDisplacementVectorsBounds2 (Mesh* mesh, Mesh* object)
{
    HEADERFUN("BuildDisplacementVectorsBounds");
    BEGIN << "Build displacement vector of boundaries : " << COLOR_BLUE << mesh->GetName () << " -> " << object->GetName () << ENDLINE;

    std::string namemesh = mesh->GetName ();
    std::string nameobject = object->GetName ();

    int numPointsMesh = mesh->GetNumberOfPoints ();
    //    int numCellsMesh = mesh->GetNumberOfCells ();
    int numEdgesMesh = mesh->GetNumberOfEdges ();
    int numCellsObject = object->GetNumberOfCells ();

    auto tagEdgeSurrogateVec = mesh->GetEdgesData ()->GetIntArrays ()->Get (nameobject + NAME_TAG_SURROGATE);

    if (tagEdgeSurrogateVec == nullptr)
    {
        ERROR << "mesh : you need to tag the edge before use tags to build a surrogate domain..." << BLINKRETURN << ENDLINE;
        return;
    }

    auto tagCellVec = object->GetCellsData ()->GetIntArrays ()->Get (namemesh + NAME_INTER);

    if (tagCellVec == nullptr)
    {
        ERROR << "object : you need to tag the cell before use tags to build a surrogate domain..." << BLINKRETURN << ENDLINE;
        return;
    }

    auto normalCell = object->GetCellsData ()->GetVecArrays ()->Get (NAME_NORMAL_ON_CELLS);

    if (normalCell == nullptr)
    {
        ERROR << "object : you need to compute normals on cells before use it..." << BLINKRETURN << ENDLINE;
        return;
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
        if (tagEdgeSurrogateVec->vec.at (std::size_t (edgeId)) == int(INTER::MIXED))
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
                    if (tagCellVec->vec.at (std::size_t (cellId)) == int(INTER::MIXED))
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

    mesh->GetPointsData ()->GetVecArrays ()->Add (nameobject + NAME_DISPLACEMENT_VECTOR, d_vec);

    ENDFUN;
    return;
}


void AddLevelSetBetween (Mesh* mesh, Mesh* object)
{
    HEADERFUN("AddLevelSetBetween");
#ifdef VERBOSE
    BEGIN << "Add level set distance : on " << COLOR_BLUE << mesh->GetName () << COLOR_DEFAULT << REVERSE << " to " << COLOR_BLUE << object->GetName () << ENDLINE;
#endif

    //    int numEdgesMesh = mesh->GetNumberOfEdges ();
    //        int numCellsMesh = mesh->GetNumberOfCells ();
    int numPointsMesh = mesh->GetNumberOfPoints ();
    //    int numPointsObject = object->GetNumberOfPoints ();
    int numCellsObject = object->GetNumberOfCells ();

    std::string nameobject = object->GetName ();
    FEStore store;
    FELocalInfos loc;
    std::vector <double> levelSetPoints (std::size_t (numPointsMesh), MAX_VALUE);
    std::vector <bool> alreadyTreated (std::size_t (numPointsMesh), false);
    int err_neg_zero = 0;
    int err_too_much = 0;

    // levelSetPoints
    for (int pointId = 0; pointId < numPointsMesh; ++pointId)
    {
        int percent = static_cast<int>(static_cast<double>(pointId + 1) / static_cast<double>(numPointsMesh) * 100.);

        COUT << "\r";
        STATUS << "compute level-set " << percent << "% ...        ";

        Point* p_cur = mesh->GetPoint (pointId);
        Cell* c_min = object->GetCell (0);
        Cell* cell = nullptr;
        double dist = MAX_VALUE;
        double distOnCell = MAX_VALUE;

        for (int idCell = 0; idCell < numCellsObject; ++idCell)
        {
            cell = object->GetCell (idCell);
            distOnCell = EuclidianDist (*cell->GetCentroid (), *p_cur);

            if (distOnCell < dist)
            {
                c_min = cell;
                dist = distOnCell;
            }
        }

        FEBase* febase = store.GetElementFor (c_min);
        febase->LocalCompute (c_min->GetCentroid (), &loc);
        Point normal = loc.normalCell;

        double value = (normal | (*p_cur - *c_min->GetCentroid ()));

        if (std::abs(value) > MAX_VALUE)
            err_too_much++;

        if (std::abs(value - NEG_ZERO) < EPSILON * EPSILON)
        {
            err_neg_zero++;
            value = 0.;
        }

        levelSetPoints.at (std::size_t (pointId)) = value;
    }

    std::cout << ENDLINE;

    mesh->GetPointsData ()->GetDoubleArrays ()->Add (nameobject + NAME_LEVELSET, levelSetPoints);

#ifdef VERBOSE
    if (err_neg_zero != 0)
        WARNING << "we found points with an anormal level-set (occure '-0') : " << err_neg_zero << ENDLINE;
    if (err_too_much != 0)
        WARNING << "we found points with an anormal level-set (occure '> " << MAX_VALUE << "') : " << err_too_much << ENDLINE;

    INFOS << "level-set computed." << ENDLINE;
    ENDFUN;
#endif

    return;
}


void TagCellsFromLevelSet (Mesh* mesh, Mesh* object)
{
    HEADERFUN ("TagCellsFromLevelSet");

#ifdef VERBOSE
    BEGIN << "tag cells from levelset entitled : " << COLOR_BLUE << object->GetName () + NAME_LEVELSET << ENDLINE;
#endif

    //    int numEdgesMesh = mesh->GetNumberOfEdges ();
    int numCellsMesh = mesh->GetNumberOfCells ();
    //    int numPointsMesh = mesh->GetNumberOfPoints ();
    //    int numPointsObject = mesh->GetNumberOfPoints ();
    std::string namemesh = mesh->GetName ();
    int count_mixed = 0;

    HetDouble::Array* levelSetPoints = mesh->GetPointsData ()->GetDoubleArrays ()->Get (object->GetName () + NAME_LEVELSET);
    //    HetInt::Array* tagintersection = mesh->GetCellsData ()->GetIntArrays ()->Get (object->GetName () + NAME_INTER);

    if (levelSetPoints == nullptr)
    {
        ERROR << "you need to compute level set on points, before use it... " << BLINKRETURN << ENDLINE;
        return;
    }

    std::vector <int> tagCellInOut (std::size_t (numCellsMesh), static_cast<int>(INTER::UNKNOWN));

    for (int i = 0; i < numCellsMesh; ++i)
    {
        Cell* c = mesh->GetCell (i);

        int* value = &tagCellInOut.at (std::size_t(i));

        int countpstv = 0;
        int countnul = 0;
        int countneg = 0;

        for (Point* p : *c->GetPoints ())
        {
            int idx = p->GetGlobalIndex ();
            double d = levelSetPoints->vec.at (std::size_t (idx));
            if (d > 1e-6)
                countpstv++;
            else if (d < -1e-6)
                countneg++;
            else
                countnul++;
        }

        if ((countneg > 0) && (countpstv > 0))
        {
            *value = static_cast<int>(INTER::MIXED);
            count_mixed++;
        }
        else if ((countneg <= 0) && (countpstv > 0))
            *value = static_cast<int>(INTER::IN);
        else
            *value = static_cast<int>(INTER::OUT);
    }

    if (count_mixed != numCellsMesh)
    {
        INFOS << mesh->GetName () << " : surrogate cell tagging is done : " << numCellsMesh << "." << ENDLINE;
    }
    else
    {
        WARNING << REVERSE << "the tagging of the surrogate cells is done but we have all the cells are mixed." << ENDLINE;

        for (int id = 0; id < numCellsMesh; ++id)
            tagCellInOut.at (std::size_t (id)) = static_cast<int>(INTER::MIXED);

    }

    mesh->GetCellsData ()->GetIntArrays ()->Add (object->GetName () + NAME_INTER, tagCellInOut);

#ifdef VERBOSE
    ENDFUN;
#endif

    return;
}

void TagEdgesFromTagCells (Mesh* mesh, Mesh* object)
{
    HEADERFUN ("TagCellsFromLevelSet");

#ifdef VERBOSE
    BEGIN << "tag edges from levelset entitled : " << COLOR_BLUE << object->GetName () + NAME_INTER << ENDLINE;
#endif

    int numEdgesMesh = mesh->GetNumberOfEdges ();
    //    int numCellsMesh = mesh->GetNumberOfCells ();
    //    int numPointsMesh = mesh->GetNumberOfPoints ();
    //    int numPointsObject = mesh->GetNumberOfPoints ();
    //    std::string nameobject = mesh->GetName ();
    int count_mixed = 0;

    auto tagCellInOut = mesh->GetCellsData ()->GetIntArrays ()->Get (object->GetName () + NAME_INTER);

    if (tagCellInOut == nullptr)
    {
        ERROR << "you need to compute tag in/out cells, before use it... " << BLINKRETURN << ENDLINE;
        return;
    }

    std::vector <int> tagEdgeSurrogate (std::size_t (numEdgesMesh), static_cast<int>(INTER::UNKNOWN));

    for (int edgeId = 0; edgeId < numEdgesMesh; ++edgeId)
    {
        Edge* edge = mesh->GetEdge (edgeId);
        int* value = &tagEdgeSurrogate.at (std::size_t (edgeId));

        std::vector <Cell*>* cellList = edge->GetCellsList ();

        int count_INSIDE = 0;
        int count_OUTSIDE = 0;
        int count_MIXED = 0;

        for (Cell* cell : *cellList)
        {
            switch (tagCellInOut->vec.at (std::size_t (cell->GetGlobalIndex ())))
            {
            case static_cast<int>(INTER::IN):
                count_INSIDE++;
                break;
            case static_cast<int>(INTER::MIXED):
                count_MIXED++;
                break;
            case static_cast<int>(INTER::OUT):
                count_OUTSIDE++;
                break;
            default:
                ERROR << "TAG UNKNOWN ??" << ENDLINE;
                break;
            }
        }

        //        ERROR << "EDGE : " << edgeId << " count cell " << cellList->size () << ENDLINE;

        if (edgeId == 34546)
        {
            //            INFOS << "count IN " << count_INSIDE << ENDLINE;
            //            INFOS << "count OUT " << count_OUTSIDE << ENDLINE;
            //            INFOS << "count MIXED " << count_MIXED << ENDLINE;
        }
        if (count_MIXED != 0 && count_OUTSIDE != 0)
        {
            *value = static_cast<int>(INTER::MIXED);
            count_mixed++;
        }
        else if ((count_INSIDE == 0) && (count_OUTSIDE != 0))
            *value = static_cast<int>(INTER::OUT);
        else
            *value = static_cast<int>(INTER::IN);
    }

    if (count_mixed != numEdgesMesh)
    {
        INFOS << mesh->GetName () << " : surrogate edge tagging is done : " << numEdgesMesh << "." << ENDLINE;
    }
    else
    {
        WARNING << REVERSE << "the tagging of the surrogate cells is done but we have all the cells are mixed." << ENDLINE;

        for (int id = 0; id < numEdgesMesh; ++id)
            tagEdgeSurrogate.at (std::size_t (id)) = static_cast<int>(INTER::MIXED);
    }

    mesh->GetEdgesData ()->GetIntArrays ()->Add (object->GetName () + NAME_TAG_SURROGATE, tagEdgeSurrogate);

#ifdef VERBOSE
    ENDFUN;
#endif

    return;
}



void BuildDisplacementVectorsBounds (Mesh* mesh, Mesh* object)
{
    HEADERFUN("BuildDisplacementVectorsBounds");
    BEGIN << "Build displacement vector of boundaries : " << COLOR_BLUE << mesh->GetName () << " -> " << object->GetName () << ENDLINE;

    std::string namemesh = mesh->GetName ();
    std::string nameobject = object->GetName ();

    int numPointsMesh = mesh->GetNumberOfPoints ();
    //    int numCellsMesh = mesh->GetNumberOfCells ();
    int numEdgesMesh = mesh->GetNumberOfEdges ();
    //    int numCellsObject = object->GetNumberOfCells ();

    // Tag Surrogate Edges Mesh (tagSEM)
    HetContainer<int>::Array* tagSEM = mesh->GetEdgesData ()->GetIntArrays ()->Get (nameobject + NAME_TAG_SURROGATE);

    if (tagSEM == nullptr)
    {
        ERROR << "mesh : you need to tag the edge before use tags to build a surrogate domain..." << BLINKRETURN << ENDLINE;
        return;
    }

    // Tag Intersection Cells Object (tagICO)
    HetContainer<int>::Array* tagICO = object->GetCellsData ()->GetIntArrays ()->Get (namemesh + NAME_INTER);

    if (tagICO == nullptr)
    {
        ERROR << "object : you need to tag the cell before use tags to build a surrogate domain..." << BLINKRETURN << ENDLINE;
        return;
    }

    // ************************************************************* //

    // Initialize new vec

    std::vector <Point*> d_vec;
    InitPointVector (&d_vec, std::size_t (numPointsMesh));

    std::vector <bool> alreadyTreated (std::size_t (numPointsMesh), false);

    // ************************************************************* //

    FELocalInfos loc;
    FEStore store;

    int count_error = 0;
    int count_d_vec = 0;

    // Loop on edges of the mesh
    for (int edgeId = 0; edgeId < numEdgesMesh; ++edgeId)
    {
        Edge* edge = mesh->GetEdge (edgeId);

        // Test if the edges is on mixed cell
        if (tagSEM->vec.at (std::size_t (edgeId)) == static_cast<int>(INTER::MIXED))
        {
            // Compute d vectors
            // if we are now, the current edge is a surrogate edge ! We need
            // to find the 'best' cell on physical domain to map all points
            // of the edge.

            for (Point* point : *edge->GetPoints ())
            {
                int pointId = point->GetGlobalIndex ();

                // this point is already treated, so skip it
                if (alreadyTreated.at (std::size_t (pointId)))
                    continue;

                Point* d = d_vec.at (std::size_t (point->GetGlobalIndex ()));
                int ret = GetDisplacementVectorAtPoint (point, object, d);

                count_d_vec++;

                if (ret == EXIT_FAILURE)
                    count_error++;

                alreadyTreated.at (std::size_t (point->GetGlobalIndex ())) = true;
            }

        }
    }

    if (count_error != 0)
        ERROR << "there is " << count_error << " \'nan\' in d vectors..." << ENDLINE;

    INFOS << "displacement vector d for points on surrogate boundary : " << count_d_vec  << ENDLINE;

    mesh->GetPointsData ()->GetVecArrays ()->Add (nameobject + NAME_DISPLACEMENT_VECTOR, d_vec);

    ENDFUN;
    return;
}


int GetDisplacementVectorAtPoint (Point* atpoint, Mesh* targetToMap, Point* d)
{
    int numCellsObject = targetToMap->GetNumberOfCells ();
    Cell* bestcell = nullptr;
    double minDist = MAX_VALUE;
    double meanDist = 0.;

    FEStore store;

    for (int cellId = 0; cellId < numCellsObject; ++cellId)
    {
        Cell* cell = targetToMap->GetCell (cellId);

        meanDist = 0.;
        for (Point* pointOfCell : *cell->GetPoints ())
            meanDist = (*pointOfCell - *atpoint).EuclidianNorm () / static_cast<double>(cell->GetPoints ()->size ());

        if (meanDist < minDist)
        {
            minDist = meanDist;
            bestcell = cell;
        }
    }

    if (bestcell == nullptr)
    {
        ERROR << "no best cell found..." << BLINKRETURN << ENDLINE;
        return EXIT_FAILURE;
    }

    FEBase* febase = store.GetElementFor (bestcell);
    *d = {0., 0., 0.};

    double coeff_norm = 1. / static_cast<double>(bestcell->GetNumberOfPoints ());

    for (Point* ptOnCell : *bestcell->GetPoints ())
    {
        Point pt = febase->TransformEleToRef (ptOnCell);
        FELocalInfos locCell;
        febase->LocalCompute (&pt, &locCell);

        Point diff = *ptOnCell - *atpoint;
        //                    *d += A - (A | loc.tangentEdge) * loc.tangentEdge;
        *d += coeff_norm * (diff | locCell.normalCell) * locCell.normalCell;
    }

    return EXIT_SUCCESS;
}
