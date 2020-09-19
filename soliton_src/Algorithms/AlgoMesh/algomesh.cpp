#include "algomesh.h"

#include <Core/core.h>
#include <IO/io.h>

void Build_NtoN (Mesh* mesh)
{
    HEADERFUN("Build_NtoN");

    //    int numpoints = mesh->GetNumberOfTotalPoints ();
    int numcells = mesh->GetNumberOfCells ();

    for (int i = 0; i < numcells; ++i)
    {
        Cell* c = mesh->GetCell (i);
        std::size_t np = std::size_t (c->GetNumberOfPoints ());
        auto lp = c->GetPoints ();

        if (c->GetTypeVTK () == VTK_CELL_TYPE::VTK_LINE ||
                c->GetTypeVTK () == VTK_CELL_TYPE::VTK_QUADRATIC_EDGE)
        {
            for (std::size_t j = 0; j < np-1; ++j)
            {
                lp->at (j)->AddPointNeighbour (lp->at (j+1));
                lp->at (j+1)->AddPointNeighbour (lp->at (j));
            }
        }
        else if (c->GetTypeVTK () == VTK_CELL_TYPE::VTK_TRIANGLE||
                 c->GetTypeVTK () == VTK_CELL_TYPE::VTK_QUADRATIC_TRIANGLE)
        {
            for (std::size_t j = 0; j < np-1; ++j)
            {
                lp->at (j)->AddPointNeighbour (lp->at (j+1));
                lp->at (j+1)->AddPointNeighbour (lp->at (j));
            }

            lp->at (np-1)->AddPointNeighbour (lp->at (0));
            lp->at (0)->AddPointNeighbour (lp->at (np-1));
        }
        else
        {
            ERROR << "the cell " << i << " (type " << c->GetTypeVTK () << ") is not supported yet.. please look the meshtools file in Build_NtoN." << ENDLINE;
        }
    }

    return;
}

void ComputeNormalsOnEdges (Mesh* mesh)
{
    HEADERFUN("ComputeNormalsOnEdges");

#ifdef VERBOSE
    BEGIN << "Compute normals on edges of the mesh : " << COLOR_BLUE << mesh->GetName () << ENDLINE;
#endif


    int numEdges = mesh->GetNumberOfEdges ();
    std::vector<Point*> vecnormals;
    vecnormals.resize (std::size_t (numEdges));

    for (int edgeId = 0; edgeId < numEdges; ++edgeId)
    {
        Edge* edge = mesh->GetEdge (edgeId);
        std::vector<Point*>* listpts = edge->GetPoints ();

        switch (edge->GetTypeVTK ())
        {
        default:
            ERROR << "only VTK_LINE and VTK_TRIANGLE are supported now to compute normals on a edge. (" << edge->GetTypeVTK () << ")" << BLINKRETURN << ENDLINE;
            return;
        case VTK_CELL_TYPE::VTK_VERTEX:
            vecnormals.at (std::size_t (edge->GetGlobalIndex ())) = new Point();
            break;
        case VTK_CELL_TYPE::VTK_LINE:
        {
            if (listpts->size () != 2)
            {
                ERROR << "something wrong, we have not 2 points for a edge, but is detected as VTK_LINE..." << BLINKRETURN << ENDLINE;
                return;
            }

            Point* normal = new Point();
            *normal = *listpts->at (0) - *listpts->at (1);
            double temp = normal->x;
            normal->x = - normal->y;
            normal->y = temp;

            *normal = *normal / normal->EuclidianNorm ();

            vecnormals.at (std::size_t (edge->GetGlobalIndex ())) = normal;

            continue;
        }
        case VTK_CELL_TYPE::VTK_TRIANGLE:
        {
            if (listpts->size () != 3)
            {
                ERROR << "something wrong, we have not 3 points for a cell, but is detected as VTK_TRIANGLE..." << BLINKRETURN << ENDLINE;
                return;
            }

            Point p1 = *listpts->at (1) - *listpts->at (0);
            Point p2 = *listpts->at (2) - *listpts->at (0);

            Point* normal = new Point();
            normal->x = p1.y * p2.z - p1.z * p2.y;
            normal->y = p1.z * p2.x - p1.x * p2.z;
            normal->z = p1.x * p2.y - p1.y * p2.x;

            *normal = *normal / normal->EuclidianNorm ();

            vecnormals.at (std::size_t (edge->GetGlobalIndex ())) = normal;

            continue;
        }
        }


    }

    mesh->GetEdgesData ()->GetVecArrays ()->Add (NAME_NORMAL_ON_EDGES, vecnormals);

#ifdef VERBOSE
    INFOS << "build : " << vecnormals.size () << " normals on edges." << ENDLINE;
    ENDFUN;
#endif

    return;
}

void ComputeNormalsOnCells (Mesh* mesh)
{
    HEADERFUN("ComputeNormalsOnCells");

#ifdef VERBOSE
    BEGIN << "Compute normals on cells of the mesh : " << COLOR_BLUE << mesh->GetName () << ENDLINE;
#endif

    int numcells = mesh->GetNumberOfCells ();
    std::vector<Point*> vecnormals;
    vecnormals.resize (std::size_t (numcells));

    for (int cellid = 0; cellid < numcells; ++cellid)
    {
        Cell* c = mesh->GetCell (cellid);
        std::vector<Point*>* listpts = c->GetPoints ();

        switch (c->GetTypeVTK ())
        {
        default:
            ERROR << "only VTK_LINE and VTK_TRIANGLE are supported now to compute normals on a cell. (" << c->GetTypeVTK () << ")" << BLINKRETURN << ENDLINE;
            return;
        case VTK_CELL_TYPE::VTK_LINE:
        {
            if (listpts->size () != 2)
            {
                ERROR << "something wrong, we have not 2 points for a cell, but is detected as VTK_LINE..." << BLINKRETURN << ENDLINE;
                return;
            }

            Point* normal = new Point();
            *normal = *listpts->at (0) - *listpts->at (1);
            double temp = normal->x;
            normal->x = - normal->y;
            normal->y = temp;

            *normal = *normal / normal->EuclidianNorm ();

            vecnormals.at (std::size_t (c->GetGlobalIndex ())) = normal;

            continue;
        }
        case VTK_CELL_TYPE::VTK_TRIANGLE:
        {
            if (listpts->size () != 3)
            {
                ERROR << "something wrong, we have not 3 points for a cell, but is detected as VTK_TRIANGLE..." << BLINKRETURN << ENDLINE;
                return;
            }

            Point p1 = *listpts->at (1) - *listpts->at (0);
            Point p2 = *listpts->at (2) - *listpts->at (0);

            Point* normal = new Point();
            normal->x = p1.y * p2.z - p1.z * p2.y;
            normal->y = p1.z * p2.x - p1.x * p2.z;
            normal->z = p1.x * p2.y - p1.y * p2.x;

            *normal = *normal / normal->EuclidianNorm ();

            vecnormals.at (std::size_t (c->GetGlobalIndex ())) = normal;

            continue;
        }
        }
    }

    mesh->GetCellsData ()->GetVecArrays ()->Add (NAME_NORMAL_ON_CELLS, vecnormals);

#ifdef VERBOSE
    INFOS << "build : " << vecnormals.size () << " normals on cells." << ENDLINE;
    ENDFUN;
#endif
    return;
}

void ComputeNormalsOnPoints(Mesh* mesh)
{
    HEADERFUN("ComputeNormalsOnPoints");

#ifdef VERBOSE
    BEGIN << "Compute normals on points of the mesh : " << COLOR_BLUE << mesh->GetName () << ENDLINE;
#endif

    int numpoints = mesh->GetNumberOfPoints ();
    std::vector<Point*> vecnormals;
    //    vecnormals.resize (std::size_t (numpoints));

    auto normalsoncells = mesh->GetCellsData ()->GetVecArrays ()->Get ("NormalOnCells");
    if (normalsoncells == nullptr)
    {
        ERROR << "you need ton compute normals on cells before to compute normals on points please ! " << BLINKRETURN << ENDLINE;
        return;
    }

    for (int ptid = 0; ptid < numpoints; ++ptid)
    {
        Point* p = mesh->GetPoint (ptid);
        std::vector<Cell*> listcells = p->GetLinkedCell ();

        Point* normal = new Point ();

        for (Cell* c : listcells)
        {
            if (c->GetCat () == CAT_CELL_EDGE::CELL)
                *normal = *normal + *normalsoncells->vec.at (std::size_t (c->GetGlobalIndex ()));
        }

        *normal = *normal / static_cast<double>(listcells.size ());

        vecnormals.push_back (normal);
    }

    mesh->GetPointsData ()->GetVecArrays ()->Add (NAME_NORMAL_ON_POINTS, vecnormals);

#ifdef VERBOSE

    INFOS << "build : " << vecnormals.size () << " normals on points." << ENDLINE;

    ENDFUN;
#endif
    return;
}

void MoveObject (Mesh* mesh, double radius, Point center)
{
    HEADERFUN ("MoveObject");
#ifdef VERBOSE
    BEGIN << "Move object " << COLOR_BLUE << mesh->GetName () << ENDLINE;
#endif

    Point mover;

    int numPoints = mesh->GetNumberOfPoints ();

    // center

    for (int pointId = 0; pointId < numPoints; ++pointId)
        mover += *mesh->GetPoint (pointId);

    mover = mover / double (numPoints);

    mover = center - mover;

#ifdef VERBOSE
    INFOS << "displacement vector    :\t" << mover << ENDLINE;
#endif

    for (int pointId = 0; pointId < numPoints; ++pointId)
    {
        Point* point = mesh->GetPoint (pointId);
        *point += mover;
    }

    // radius
    mover = center + Point({0, 0, mover.z});

#ifdef VERBOSE
    INFOS << "new center             :\t" << mover << ENDLINE;
#endif

    double maxradius = 0;
    for (int pointId = 0; pointId < numPoints; ++pointId)
    {
        double radiustemp = EuclidianDist (mover, *mesh->GetPoint (pointId));

        if (radiustemp > maxradius)
            maxradius = radiustemp;
    }

#ifdef VERBOSE
    INFOS << "current maximum radius :\t" << maxradius << ENDLINE;
    INFOS << "desired radius         :\t" << radius << ENDLINE;
#endif

    double coeff  = radius / maxradius;

#ifdef VERBOSE
    INFOS << "scaling coefficient    :\t" << coeff << ENDLINE;
#endif

    for (int pointId = 0; pointId < numPoints; ++pointId)
    {
        Point* point = mesh->GetPoint (pointId);
        Point distorcer = *point - mover;

        distorcer *= coeff;
        *point = mover + distorcer;
    }

#ifdef VERBOSE
    STATUS << "displacement done !" << ENDLINE;
    ENDFUN;
#endif

    return;
}

void TransfertTagPhysical (Mesh* mesh)
{
    HEADERFUN ("TransfertTagPhysical");

    auto tagvec = mesh->GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec == nullptr)
    {
        INFOS << "no tag physical on cells : search under the name " << NAME_TAG_PHYSICAL << "! " << BLINKRETURN << ENDLINE;
        return;
    }

    int numPoints = mesh->GetNumberOfPoints ();
    int numCells = mesh->GetNumberOfCells ();
    int numEdges = mesh->GetNumberOfEdges ();

    std::vector<int> tagOnPt (std::size_t (numPoints), static_cast<int>(TAG_PHYSICAL::TAG_DOMAIN));

    for (int i = 0; i < numCells; ++i)
    {
        TAG_PHYSICAL tag = static_cast<TAG_PHYSICAL>(tagvec->vec.at (std::size_t (i)));

        if (tag == TAG_PHYSICAL::TAG_DOMAIN)
            continue;

        for (Point* p : *mesh->GetCell (i)->GetPoints ())
            tagOnPt.at (std::size_t (p->GetGlobalIndex ())) = static_cast<int>(tag);
    }

    mesh->GetPointsData ()->GetIntArrays ()->Add (NAME_TAG_PHYSICAL, tagOnPt);

    std::vector<int> tagOnEdges (std::size_t (numEdges), static_cast<int>(TAG_PHYSICAL::TAG_DOMAIN));

    for (int i = 0; i < numEdges; ++i)
    {
        std::vector<Point*>* listPts = mesh->GetEdge (i)->GetPoints ();

        int count_wall = 0;
        int count_inlet = 0;
        int count_outlet = 0;
        int count_domain = 0;

        for (Point* p : *listPts)
        {
            switch (TAG_PHYSICAL (tagOnPt.at (std::size_t (p->GetGlobalIndex ()))))
            {
            case TAG_PHYSICAL::TAG_WALL:
                count_wall++;
                break;
            case TAG_PHYSICAL::TAG_INLET:
                count_inlet++;
                break;
            case TAG_PHYSICAL::TAG_OUTLET:
                count_outlet++;
                break;
            case TAG_PHYSICAL::TAG_DOMAIN:
                count_domain++;
                break;
            default:
                break;
            }
        }

        if (count_domain != 0)
            continue;

        TAG_PHYSICAL tag = TAG_PHYSICAL::TAG_DOMAIN;

        if (count_wall != 0)
            tag = TAG_PHYSICAL::TAG_WALL;
        if (count_inlet != 0)
            tag = TAG_PHYSICAL::TAG_INLET;
        if (count_outlet != 0)
            tag = TAG_PHYSICAL::TAG_OUTLET;


        tagOnEdges.at (std::size_t (i)) = static_cast<int>(tag);
    }

    mesh->GetEdgesData ()->GetIntArrays ()->Add (NAME_TAG_PHYSICAL, tagOnEdges);

    return;
}


void ComputeDampingArea (Mesh* mesh, TAG_PHYSICAL tag, double h)
{
    int numPoints = mesh->GetNumberOfPoints ();
    int numCells = mesh->GetNumberOfCells ();

    std::vector<bool> tagpoint (std::size_t(numPoints), false);
    std::vector<bool> tagcell (std::size_t(numCells), false);


    for (int pointId = 0; pointId < numPoints; ++pointId)
    {
        double value = GetDampingCoeffFor (mesh->GetPoint (pointId), mesh, tag, h);

        if (value >= 0)
            tagpoint.at (std::size_t (pointId)) = true;
    }

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        bool value = true;
        Cell* cell = mesh->GetCell (cellId);

        for (Point* pointOnCell : *cell->GetPoints ())
        {
            int pointId = pointOnCell->GetGlobalIndex ();

            if (!tagpoint.at (std::size_t (pointId)))
            {
                value = false;
                break;
            }
        }

        if (value)
            tagcell.at (std::size_t (cellId)) = true;
    }

    mesh->GetPointsData ()->GetBooleanArrays ()->Add (NAME_TAG_DAMPING_AREA, tagpoint);
    mesh->GetCellsData ()->GetBooleanArrays ()->Add (NAME_TAG_DAMPING_AREA, tagcell);

    return;
}

double GetDampingCoeffFor (Point* atpoint, Mesh* mesh, TAG_PHYSICAL tag, double h)
{
    int numPoints = mesh->GetNumberOfPoints ();
    HetInt::Array* tagphysical = mesh->GetPointsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    double distmin = 1e6;

    for (int pointId = 0; pointId < numPoints; ++pointId)
    {
        if (tagphysical->vec.at (std::size_t (pointId)) != static_cast<int>(tag))
            continue;
        Point* point = mesh->GetPoint (pointId);

        double dist = EuclidianDist (*point, *atpoint);

        if (dist < distmin)
            distmin = dist;
    }


    if (distmin < std::abs(10. * h))
        return distmin;
    return -1.;

}
