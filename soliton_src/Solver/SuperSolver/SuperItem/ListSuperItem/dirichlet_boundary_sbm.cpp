#include <Solver/SuperSolver/SuperItem/superitem.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>
#include <Algorithms/algorithms.h>

/** ********************************************************************** *
 * ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM
 * *********************************************************************** *
 *
 * d    = displacement vector
 * [C]  = coeff penalization
 * U_D  = U on boundary
 * {X}  = surrogate resp of X
 *
 *
 * MATRIX :
 * -1 <phi_j | grad_phi_i.{normal}>
 * -1 <grad_phi_j.{normal} | phi_i + grad_phi_i.d>
 * +1 <[C] * (phi_j + grad_phi_j.d) | phi_i + grad_phi_i.d>
 *
 * SECOND MEMBRE :
 * -1 <grad_phi_j.{normal} | {U_D}>
 * +1 <[C] * (phi_j + grad_phi_j.d) | {U_D}>
 **/

template<>
void
SuperItem::InternalCompute<ITEM_T::DIRICHLET_BOUNDARY_SBM> ()
{
    SuperItem *itemThis = this;

    if (!(itemThis->ItsVariableOverTime ()|| itemThis->GetForce ()))
        return;

    std::vector<Duo> duos;
    std::vector<Triplet> triplets;
    Mesh* mesh = itemThis->GetMesh ();
    int numCells  = itemThis->GetMesh ()->GetNumberOfCells ();
    int szDisplayNumCells = static_cast<int>(std::to_string (numCells).size ());
    int numPoints = itemThis->GetMesh ()->GetNumberOfPoints ();
    int numEdges  = itemThis->GetMesh ()->GetNumberOfEdges ();

    FEStore *festore = new FEStore ();
    QuadStore *quadstore = new QuadStore ();
    duos.reserve (static_cast<std::size_t>(32 * numPoints));
    triplets.reserve (static_cast<std::size_t>(32 * numPoints));

    if (itemThis->GetCollect ())
    {
        itemThis->GetApplicationCell ()->resize (static_cast<std::size_t>(numCells), 0);
        itemThis->GetApplicationEdge ()->resize (static_cast<std::size_t>(numEdges), 0);
        itemThis->GetApplicationPoint ()->resize (static_cast<std::size_t>(numPoints), 0);
    }

    //    HetContainer<int>::Array* listTagsPhysicalCell = mesh->GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);
    HetContainer<int>::Array* listTagsPhysicalEdge = mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);
    HetContainer<int>::Array* listTagsSurrogateCell = mesh->GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_INTERSECTION);
    HetContainer<int>::Array* listTagsSurrogateEdge = mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_INTERSECTION);

    HetContainer<Point*>::Array* listDVec = mesh->GetPointsData ()->GetVecArrays ()->Get (NAME_DISPLACEMENT_VECTOR);


    double penal = itemThis->GetCoeffPen () / (std::pow(mesh->GetPrescribedSize (), itemThis->GetPowPen ()));

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell                    *cell           = mesh->GetCell (cellId); \
        FEBase                  *febase4Cell    = festore->GetElementFor (cell);
        //        QuadStore::QuadObject   *quad4Cell      = quadstore->Get (febase4Cell);
        std::size_t             numEdgesOnCell  = std::size_t (cell->GetNumberOfEdges ());
        std::size_t             numPointsOnCell = std::size_t (cell->GetNumberOfPoints());

        if (itemThis->GetView ())
        {
            int percent = static_cast<int>(static_cast<double>(cellId + 1) / static_cast<double>(numCells) * 100.); \
            std::cout << "\r\t";  \
            TREE_BRANCH << COLOR_BLUE       << "[" << std::setw(3) << percent << "%] " << std::flush;
            COUT        << COLOR_DEFAULT    << *itemThis << " " << std::flush;
            COUT        << COLOR_MAGENTA    << "[Cell " << std::flush;
            COUT        << COLOR_BLUE       << std::setw(szDisplayNumCells) << cellId + 1 << "/" << numCells << std::flush;
            COUT        << COLOR_MAGENTA    << " " << to_string(cell->GetTypeVTK()) << "]                         " << std::flush;
        }

        INTER tagSurrogateCell = static_cast<INTER>(listTagsSurrogateCell->vec.at (static_cast<std::size_t>(cell->GetGlobalIndex ())));
        if (tagSurrogateCell != itemThis->GetSolver ()->GetInterTag () &&
                tagSurrogateCell != INTER::DEFAULT &&
                itemThis->GetSolver ()->GetInterTag () != INTER::DEFAULT)
            continue;

        //        int tagPhysicalCell = listTagsPhysicalCell->vec.at (static_cast<std::size_t>(cell->GetGlobalIndex ()));

        //        if (tagPhysicalCell != static_cast<int>(itemThis->GetTag2Apply ()) && tagPhysicalCell != static_cast<int>(PHYS::DEFAULT))
        //            continue;

        if (itemThis->GetCollect ())\
            itemThis->GetApplicationCell ()->at (static_cast<std::size_t>(cellId)) += 1;

        for (std::size_t edgeId = 0; edgeId < numEdgesOnCell; ++edgeId)
        {
            Edge*                   edge            = cell->GetEdges()->at(edgeId);
            Edge*                   edgeOnCellRef   = febase4Cell->GetEdge (edgeId);
            FEBase*                 febase4Edge     = festore->GetElementFor (edgeOnCellRef);
            QuadStore::QuadObject*  quad4Edge       = quadstore->Get (febase4Edge);


            int tagPhysicalEdge= listTagsPhysicalEdge->vec.at (static_cast<std::size_t>(edge->GetGlobalIndex ()));

            if (tagPhysicalEdge!= static_cast<int>(itemThis->GetTag2Apply ()) &&
                    tagPhysicalEdge != static_cast<int>(PHYS::DEFAULT) &&
                    static_cast<int>(itemThis->GetTag2Apply ()) != static_cast<int>(PHYS::DEFAULT))
                continue;


            int tagSurrogateEdge = listTagsSurrogateEdge->vec.at (static_cast<std::size_t>(edge->GetGlobalIndex ()));

            if (tagSurrogateEdge != static_cast<int>(INTER::MIXED))
                continue;

            if (itemThis->GetCollect ())
                itemThis->GetApplicationEdge ()->at (static_cast<std::size_t>(edge->GetGlobalIndex ())) += 1;

            for (std::size_t i = 0; i < static_cast<std::size_t>(numPointsOnCell); ++i)
            {
                LambdaOnPoint<double>   fun_phi_i       = febase4Cell->GetPhi (i);
                LambdaOnPoint<Point>    fun_grad_phi_i  = febase4Cell->GetGradPhi (i);
                int                     id_i            = cell->GetPoints ()->at (i)->GetGlobalIndex ();

                if (itemThis->GetCollect ())
                    itemThis->GetApplicationPoint ()->at (static_cast<std::size_t>(id_i)) += 1;

                for (std::size_t j = 0; j < static_cast<std::size_t>(numPointsOnCell); ++j)
                {
                    LambdaOnPoint<double>   fun_phi_j       = febase4Cell->GetPhi (j);
                    LambdaOnPoint<Point>    fun_grad_phi_j  = febase4Cell->GetGradPhi (j);
                    int                     id_j            = cell->GetPoints ()->at (j)->GetGlobalIndex ();

                    if (itemThis->GetCollect ())
                        itemThis->GetApplicationPoint ()->at (static_cast<std::size_t>(id_j)) += 1;


                    double coeff_matrix = 0.;
                    double coeff_sec    = 0.;
                    Point dvec;

                    for (Point *ptOnEdge : *edge->GetPoints ())
                        dvec += *listDVec->vec.at (static_cast<std::size_t>(ptOnEdge->GetGlobalIndex ()));

                    dvec /= static_cast<double>(edge->GetNumberOfPoints ());

//                    *listOfUsefulDVec.at (static_cast<std::size_t>(edge->GetGlobalIndex ())) = dvec;

                    for (std::size_t k = 0; k < quad4Edge->npts; ++k)
                    {
                        Point point_int     = quad4Edge->pts [k];
                        Point point_tran    = febase4Edge->TransformRefToEle (&point_int);
                        Point point_real    = febase4Cell->TransformRefToEle (&point_tran);

                        FELocalInfos locEdge; febase4Edge->LocalCompute (&point_int, &locEdge);
                        FELocalInfos locCell; febase4Cell->LocalCompute (&point_tran, &locCell);

                        double weight = locEdge.detJac * locCell.detJac * quad4Edge->w [k];

                        Point grad_phi_i = locEdge.JacInvT * locCell.JacInvT * fun_grad_phi_i (&point_tran);
                        Point grad_phi_j = locEdge.JacInvT * locCell.JacInvT * fun_grad_phi_j (&point_tran);

                        double phi_i = fun_phi_i (&point_tran);
                        double phi_j = fun_phi_j (&point_tran);

                        double fun   = itemThis->GetFun () (point_real, itemThis->GetTime ());
                        Point normal = locEdge.normalCell_ref;


                        Point d = dvec;
                        d = locCell.Jac_cmp * locEdge.Jac_cmp * d;

                        //!! -1 <phi_j | grad_phi_i.{normal}>
                        coeff_matrix -= weight * phi_j * (grad_phi_i | normal);

                        //!! -1 <grad_phi_j.{normal} | phi_i + grad_phi_i.d>
                        coeff_matrix -= weight * (grad_phi_j | normal) * (phi_i + (grad_phi_i | d));

                        //!! +1 <[C] * (phi_j + grad_phi_j.d) | phi_i + grad_phi_i.d>
                        coeff_matrix += weight * penal * (phi_j + (grad_phi_j | d)) * (phi_i + (grad_phi_i | d));

                        if (i == 0)
                        {
                            //!! -1 <grad_phi_j.{normal} | {U_D}>
                            coeff_sec -= weight * (grad_phi_j | normal) * fun;

                            //!! +1 <[C] * (phi_j + grad_phi_j.d) | {U_D}>
                            coeff_sec += weight * (phi_j + (grad_phi_j | d)) * fun;
                        }

                    }

                    triplets.push_back (Triplet (id_i, id_j, coeff_matrix));
                    if (i == 0)
                        duos.push_back (Duo (id_j, coeff_sec));
                }
            }
        }
    }

//    mesh->GetEdgesData ()->GetVecArrays ()->Add ("d_vec" +  itemThis->GetWholeName (), listOfUsefulDVec);

    itemThis->SetSparseMatrix (&triplets);
    itemThis->SetSecMember (&duos);

    delete festore;
    delete quadstore;

    if (itemThis->GetCollect ())
    {
        HetContainer<int>::Array* arrayCell = mesh->GetCellsData ()->GetIntArrays ()->Get (itemThis->GetWholeName () +  NAME_CONTRIBUTIONS_ON_CELLS);
        if (arrayCell != nullptr)
            arrayCell->vec = *itemThis->GetApplicationCell ();
        else
            mesh->GetCellsData ()->GetIntArrays ()->Add (itemThis->GetWholeName () +  NAME_CONTRIBUTIONS_ON_CELLS, *itemThis->GetApplicationCell ());

        HetContainer<int>::Array* arrayEdge = mesh->GetEdgesData ()->GetIntArrays ()->Get (itemThis->GetWholeName () +  NAME_CONTRIBUTIONS_ON_EDGES);
        if (arrayEdge != nullptr)
            arrayEdge->vec = *itemThis->GetApplicationEdge ();
        else
            mesh->GetEdgesData ()->GetIntArrays ()->Add (itemThis->GetWholeName () +  NAME_CONTRIBUTIONS_ON_EDGES, *itemThis->GetApplicationEdge ());

        HetContainer<int>::Array* arrayPoint = mesh->GetPointsData ()->GetIntArrays ()->Get (itemThis->GetWholeName () +  NAME_CONTRIBUTIONS_ON_POINTS);
        if (arrayPoint != nullptr)
            arrayPoint->vec = *itemThis->GetApplicationPoint ();
        else
            mesh->GetPointsData ()->GetIntArrays ()->Add (itemThis->GetWholeName () +  NAME_CONTRIBUTIONS_ON_POINTS, *itemThis->GetApplicationPoint ());
    }

    if (itemThis->GetView ())
        ENDFUN;

    return;
}