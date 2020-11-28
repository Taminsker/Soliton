#include <Solver/SolitonFE/SolitonFEItem/solitonfeitem.h>
#include <Solver/SolitonFE/SolitonFEContainer/solitonfecontainer.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>
#include <Algorithms/algorithms.h>

/** ********************************************************************** *
 * ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY
 * ********************************************************************** *
 *
 * U_D = U on boundary
 * [C] = coeff penalization
 *
 * MATRIX :
 * -1 <phi_j, grad_phi_i.normal>
 * -1 <grad_phi_j.normal | phi_i>
 * +1 <[C] * phi_j | phi_i>
 *
 * SECOND MEMBRE :
 * -1 <grad_phi_j.normal | U_D>
 * +1 <[C] * phi_j | U_D>
 **/

template<>
void
SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY>::Compute (real_t time)
{
    if (!(m_var0vrT || m_c->m_force_compute))
        return;

    std::vector<Duo> duos;
    std::vector<Triplet_eig> triplets;
    Mesh* mesh = m_c->m_store->GetMainMesh ();
    int numCells = mesh->GetNumberOfCells ();
    int szDisplayNumCells = static_cast<int>(std::to_string (numCells).size ());
    int numPoints = mesh->GetNumberOfPoints ();
    int numEdges = mesh->GetNumberOfEdges ();

    FEStore *festore = new FEStore ();
    QuadStore *quadstore = new QuadStore ();
    duos.reserve (static_cast<ul_t>(32 * numPoints));
    triplets.reserve (static_cast<ul_t>(32 * numPoints));

    if (m_c->m_collect_contributions)
    {
        m_appCells->resize (static_cast<ul_t>(numCells), 0);
        m_appEdges->resize (static_cast<ul_t>(numEdges), 0);
        m_appPoints->resize (static_cast<ul_t>(numPoints), 0);
    }

    //  HetContainer<int>::Array* listTagsPhysicalCell = mesh->GetCellsData ()->Get<int> (NAME_TAG_PHYSICAL);
    HetContainer<int>::Array* listTagsPhysicalEdge = mesh->GetEdgesData ()->Get<int> (NAME_TAG_PHYSICAL);
    HetContainer<int>::Array* listTagsSurrogateCell = mesh->GetCellsData ()->Get<int> (NAME_TAG_INTERSECTION);
    //  HetContainer<int>::Array* listTagsSurrogateEdge = mesh->GetEdgesData ()->Get<int> (NAME_TAG_INTERSECTION);

    //  if (listTagsPhysicalCell == nullptr || listTagsPhysicalEdge == nullptr)
    //  {
    //    ERROR << "tags physical cell is at the ptr " << listTagsPhysicalCell << " item " << to_string(nameTypeItem) << ENDLINE;
    //    ERROR << "tags physical edge is at the ptr " << listTagsPhysicalEdge << " item " << to_string(nameTypeItem) << ENDLINE;

    //    return;
    //  }

    real_t penal = m_c->m_coeff_pen / (std::pow(mesh->GetPrescribedSize (), m_c->m_pow_pen));

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell          *cell      = mesh->GetCell (cellId); \
        FEBase         *febase4Cell  = festore->GetElementFor (cell);
        //    QuadStore::QuadObject  *quad4Cell   = quadstore->Get (febase4Cell);
        ul_t       numEdgesOnCell = ul_t (cell->GetNumberOfEdges ());
        ul_t       numPointsOnCell = ul_t (cell->GetNumberOfPoints());

        if (m_c->m_view)
        {
            int percent = static_cast<int>(static_cast<real_t>(cellId + 1) / static_cast<real_t>(numCells) * 100.); \
            std::cout << "\r\t"; \
            TREE_BRANCH << COLOR_BLUE    << "[" << std::setw(3) << percent << "%] " << std::flush;
            this->Print (COUT);
            COUT    << COLOR_MAGENTA  << "[Cell " << std::flush;
            COUT    << COLOR_BLUE    << std::setw(szDisplayNumCells) << cellId + 1 << "/" << numCells << std::flush;
            COUT    << COLOR_MAGENTA  << " " << to_string(cell->GetTypeVTK()) << "]             " << std::flush;
        }

        int tagSurrogateCell = listTagsSurrogateCell->vec.at (static_cast<ul_t>(cellId));

        if (tagSurrogateCell != int(m_c->GetInterTag ()) && m_c->GetInterTag () != INTER::DEFAULT)
            continue;

        //    int tagPhysicalCell = listTagsPhysicalCell->vec.at (static_cast<ul_t>(cell->GetGlobalIndex ()));

        //    if (tagPhysicalCell != static_cast<int>(m_tag2Apply) && tagPhysicalCell != static_cast<int>(PHYS::DEFAULT))
        //      continue;

        if (m_c->m_collect_contributions)
            m_appCells->at (static_cast<ul_t>(cellId)) += 1;

        for (ul_t edgeId = 0; edgeId < numEdgesOnCell; ++edgeId)
        {
            Edge*          edge      = cell->GetEdges()->at(edgeId);
            Edge*          edgeOnCellRef  = febase4Cell->GetEdge (edgeId);
            FEBase*         febase4Edge   = festore->GetElementFor (edgeOnCellRef);
            QuadStore::QuadObject* quad4Edge    = quadstore->Get (febase4Edge);


            int tagPhysicalEdge= listTagsPhysicalEdge->vec.at (static_cast<ul_t>(edge->GetGlobalIndex ()));

            if (tagPhysicalEdge!= static_cast<int>(m_tag2Apply) && tagPhysicalEdge != static_cast<int>(PHYS::DEFAULT))
                continue;

            if (m_c->m_collect_contributions)
                m_appEdges->at (static_cast<ul_t>(edge->GetGlobalIndex ())) += 1;

            for (ul_t i = 0; i < static_cast<ul_t>(numPointsOnCell); ++i)
            {
                LambdaOnPoint<real_t>  fun_phi_i    = febase4Cell->GetPhi (i);
                LambdaOnPoint<Point>  fun_grad_phi_i = febase4Cell->GetGradPhi (i);
                int           id_i      = cell->GetPoints ()->at (i)->GetGlobalIndex ();

                if (m_c->m_collect_contributions)
                    m_appPoints->at (static_cast<ul_t>(id_i)) += 1;

                for (ul_t j = 0; j < static_cast<ul_t>(numPointsOnCell); ++j)
                {
                    LambdaOnPoint<real_t>  fun_phi_j    = febase4Cell->GetPhi (j);
                    LambdaOnPoint<Point>  fun_grad_phi_j = febase4Cell->GetGradPhi (j);
                    int           id_j      = cell->GetPoints ()->at (j)->GetGlobalIndex ();

                    if (m_c->m_collect_contributions)
                        m_appPoints->at (static_cast<ul_t>(id_j)) += 1;

                    real_t coeff_matrix = 0.;
                    real_t coeff_sec  = 0.;

                    for (ul_t k = 0; k < quad4Edge->npts; ++k)
                    {
                        Point point_int   = quad4Edge->pts [k];
                        Point point_tran  = febase4Edge->TransformRefToEle (&point_int);
                        Point point_real  = febase4Cell->TransformRefToEle (&point_tran);

                        FELocalInfos locEdge; febase4Edge->LocalCompute (&point_int, &locEdge);
                        FELocalInfos locCell; febase4Cell->LocalCompute (&point_tran, &locCell);

                        real_t weight = locEdge.detJac * locCell.detJac * quad4Edge->w [k];

                        Point grad_phi_i = locEdge.JacInvT * locCell.JacInvT * fun_grad_phi_i (&point_tran);
                        Point grad_phi_j = locEdge.JacInvT * locCell.JacInvT * fun_grad_phi_j (&point_tran);

                        real_t phi_i = fun_phi_i (&point_tran);
                        real_t phi_j = fun_phi_j (&point_tran);

                        real_t fun  = m_fun->ToReal (&point_real, time, cell);
                        Point normal = locEdge.normalCell_ref;


                        //!! -1 <phi_j, grad_phi_i.normal>
                        coeff_matrix -= weight * phi_j * (grad_phi_i | normal);
                        //!! -1 <grad_phi_j.normal | phi_i>
                        coeff_matrix -= weight * (grad_phi_j | normal) * phi_i;
                        //!! +1 <[C] * phi_j | phi_i>
                        coeff_matrix += weight * penal * phi_j * phi_i;

                        if (i == 0)
                        {
                            //!! -1 <grad_phi_j.normal | U_D>
                            coeff_sec -= weight * (grad_phi_j | normal) * fun;
                            //!! +1 <[C] * phi_j | U_D>
                            coeff_sec += weight * penal * phi_j * fun;
                        }

                    }

                    triplets.push_back (Triplet_eig (id_i, id_j, coeff_matrix));
                    if (i == 0)
                        duos.push_back (Duo (id_j, coeff_sec));
                }
            }
        }
    }

    this->SetSparseMatrix (&triplets);
    this->SetSecMember (&duos);

    delete festore;
    delete quadstore;


    if (m_c->m_collect_contributions)
    {
        HetContainer<int>::Array* arrayCell = mesh->GetCellsData ()->Get<int> (m_name + NAME_CONTRIBUTIONS_ON_CELLS);
        if (arrayCell != nullptr)
            arrayCell->vec = *m_appCells;
        else
            mesh->GetCellsData ()->Add<int> (m_name + NAME_CONTRIBUTIONS_ON_CELLS, *m_appCells);

        HetContainer<int>::Array* arrayEdge = mesh->GetEdgesData ()->Get<int> (m_name + NAME_CONTRIBUTIONS_ON_EDGES);
        if (arrayEdge != nullptr)
            arrayEdge->vec = *m_appEdges;
        else
            mesh->GetEdgesData ()->Add<int> (m_name + NAME_CONTRIBUTIONS_ON_EDGES, *m_appEdges);

        HetContainer<int>::Array* arrayPoint = mesh->GetPointsData ()->Get<int> (m_name + NAME_CONTRIBUTIONS_ON_POINTS);
        if (arrayPoint != nullptr)
            arrayPoint->vec = *m_appPoints;
        else
            mesh->GetPointsData ()->Add<int> (m_name + NAME_CONTRIBUTIONS_ON_POINTS, *m_appPoints);
    }

    if (m_c->m_view)
        ENDFUN;

    return;
}
