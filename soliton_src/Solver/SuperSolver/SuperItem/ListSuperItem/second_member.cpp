#include <Solver/SuperSolver/SuperItem/superitem.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>
#include <Algorithms/algorithms.h>


/** ********************************************************************** *
  * ITEM_SOLVER_TYPE::SECOND_MEMBER
  * ********************************************************************** *
  *  +<F, PHI(J)>
  **/
template<>
void
SuperItem::InternalCompute<ITEM_T::SECOND_MEMBER> ()
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

    HetContainer<int>::Array* listTagsPhysicalCell = mesh->GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);
//    HetContainer<int>::Array* listTagsPhysicalEdge = mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);
    HetContainer<int>::Array* listTagsSurrogateCell = mesh->GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_INTERSECTION);
//    HetContainer<int>::Array* listTagsSurrogateEdge = mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_INTERSECTION);

//    if (listTagsPhysicalCell == nullptr || listTagsPhysicalEdge == nullptr)
//    {
//        ERROR << "tags physical cell is at the ptr " << listTagsPhysicalCell << " item " << to_string(nameTypeItem) << ENDLINE;
//        ERROR << "tags physical edge is at the ptr " << listTagsPhysicalEdge << " item " << to_string(nameTypeItem) << ENDLINE;

//        return;
//    }

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell                    *cell           = mesh->GetCell (cellId); \
        FEBase                  *febase4Cell    = festore->GetElementFor (cell);
        QuadStore::QuadObject   *quad4Cell      = quadstore->Get (febase4Cell);
        //        std::size_t             numEdgesOnCell  = std::size_t (cell->GetNumberOfEdges ());
        std::size_t             numPointsOnCell = std::size_t (cell->GetNumberOfPoints());

        if (itemThis->GetView ())
        {
            int percent = static_cast<int>(static_cast<double>(cellId + 1) / static_cast<double>(numCells) * 100.); \
            std::cout << "\r\t";  \
            TREE_BRANCH << COLOR_BLUE       << "[" << std::setw(3) << percent << "%] " << std::flush;
            COUT        << COLOR_DEFAULT    << *itemThis << " " << std::flush;
            COUT        << COLOR_MAGENTA    << "[Cell " << std::flush;
            COUT        << COLOR_BLUE       << std::setw (szDisplayNumCells) << cellId + 1 << "/" << numCells << std::flush;
            COUT        << COLOR_MAGENTA    << " " << to_string(cell->GetTypeVTK()) << "]                         " << std::flush;
        }

        INTER tagSurrogateCell = static_cast<INTER>(listTagsSurrogateCell->vec.at (static_cast<std::size_t>(cell->GetGlobalIndex ())));
        if (tagSurrogateCell != itemThis->GetSolver ()->GetInterTag () &&
                tagSurrogateCell != INTER::DEFAULT &&
                itemThis->GetSolver ()->GetInterTag () != INTER::DEFAULT)
            continue;

        int tagPhysicalCell = listTagsPhysicalCell->vec.at (static_cast<std::size_t>(cell->GetGlobalIndex ()));

        if (tagPhysicalCell != static_cast<int>(itemThis->GetTag2Apply ()) && tagPhysicalCell != static_cast<int>(PHYS::DEFAULT))
            continue;

        if (itemThis->GetCollect ())\
            itemThis->GetApplicationCell ()->at (static_cast<std::size_t>(cellId)) += 1;

        std::size_t i = 0;
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

//            double coeff_matrix = 0.;
            double coeff_sec    = 0.;

            for (std::size_t k = 0; k < quad4Cell->npts; ++k)
            {
                Point        point_int  = quad4Cell->pts [k];
                Point        point_real = febase4Cell->TransformRefToEle (&point_int);
                FELocalInfos locCell; febase4Cell->LocalCompute (&point_int, &locCell);

                double weight = locCell.detJac * quad4Cell->w [k];

//                Point grad_phi_i = locCell.JacInvT * fun_grad_phi_i (&point_int);
//                Point grad_phi_j = locCell.JacInvT * fun_grad_phi_j (&point_int);

//                double phi_i = fun_phi_i (&point_int);
                double phi_j = fun_phi_j (&point_int);

                double fun = itemThis->GetFun() (point_real, itemThis->GetTime ());

                coeff_sec += weight * fun * phi_j;
            }


            duos.push_back (Duo (id_j, coeff_sec));
        }

    }

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
