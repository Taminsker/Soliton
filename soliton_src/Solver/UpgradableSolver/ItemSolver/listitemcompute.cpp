//#include <Solver/UpgradableSolver/ItemSolver/itemsolver.h>
//#include <Solver/FEStruct/festruct.h>
//#include <Solver/QuadStruct/quadstruct.h>
//#include <Algorithms/algorithms.h>

//#include "def_macros_items.h"

//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::EMPTY>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    return;
//}

///** ********************************************************************** *
//  * ITEM_SOLVER_TYPE::PHI_PHI
//  * ********************************************************************** *
//  *
//  * +<PHI(J), PHI(I)>
//  *
//  **/
//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::PHI_PHI>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}

//        if (tagPhysicalCell != m_tagToApply)
//            continue;

//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                double COEFF_MATRIX = 0.;
//                double COEFF_SEC    = 0.;

//                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
//                {
//                    OBJECTS_ON_QUAD_CELL(K, I, J);

//                    COEFF_MATRIX += WEIGHT * PHI(J) * PHI(I);
//                    if (I == 0)
//                        COEFF_SEC += WEIGHT * 0.;
//                }

//                m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});
//                if (I == 0)
//                    m_duos.push_back ({ID(J), COEFF_SEC});
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    return;
//}

///** ********************************************************************** *
//  * ITEM_SOLVER_TYPE::DAMPING_PHI_PHI
//  * ********************************************************************** *
//  *  +<PHI(J), chi PHI(I)>
//  **/
//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::DAMPING_PHI_PHI>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    auto FunDamp = [this](Point* p) -> double
//    {return GetDampingCoeffFor (p, *m_mesh, m_tagToApply, (*m_coeff_pen));};

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}


//        if (tagPhysicalCell != m_tagToApply)
//            continue;

//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                double COEFF_MATRIX = 0.;
//                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
//                {
//                    OBJECTS_ON_QUAD_CELL(K, I, J);

//                    COEFF_MATRIX += WEIGHT *  1. / FunDamp (&POINT_REAL) * PHI(J) * PHI(I);
//                }

//                m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    return;
//}

///** ********************************************************************** *
//  * ITEM_SOLVER_TYPE::GRAD_GRAD
//  * ********************************************************************** *
//  *  +<GRAD_PHI(J), GRAD_PHI(I)>
//  **/
//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::GRAD_GRAD>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();
//    int numEdges  = (*m_mesh)->GetNumberOfEdges ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    if (*m_collect)
//    {
//        m_applicableCells.clear ();
//        m_applicableEdges.clear ();

//        m_applicableCells.resize (static_cast<std::size_t>(numCells), 0x0);
//        m_applicableEdges.resize (static_cast<std::size_t>(numEdges), 0x0);
//    }

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}

//        if (tagPhysicalCell != m_tagToApply)
//            continue;

//        if (*m_collect)
//            m_applicableCells.at (static_cast<std::size_t>(CELLID)) += 1;

//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                double COEFF_MATRIX = 0.;

//                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
//                {
//                    OBJECTS_ON_QUAD_CELL(K, I, J);

//                    COEFF_MATRIX += WEIGHT * (GRAD_PHI(J) | GRAD_PHI(I));
//                }

//                m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    return;
//}

///** ********************************************************************** *
//  * ITEM_SOLVER_TYPE::SECOND_MEMBER
//  * ********************************************************************** *
//  *  +<F, PHI(J)>
//  **/
//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::SECOND_MEMBER>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    if (*m_collect)
//    {
//        m_applicableCells.clear ();
//        m_applicableEdges.clear ();
//    }

//    int numCells  = (*m_mesh)->GetNumberOfCells ();
//    int numEdges  = (*m_mesh)->GetNumberOfEdges ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    if (*m_collect)
//    {
//        m_applicableCells.resize (static_cast<std::size_t>(numCells), 0);
//        m_applicableEdges.resize (static_cast<std::size_t>(numEdges), 0);
//    }

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}

//        if (tagPhysicalCell != m_tagToApply)
//            continue;

//        if (*m_collect)
//            m_applicableCells.at (static_cast<std::size_t>(CELLID)) += 1;

//        std::size_t I = 0;
//        OBJECTS_ON_POINT(I);

//        for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//        {
//            OBJECTS_ON_POINT(J);

//            double COEFF_SEC    = 0.;

//            for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
//            {
//                OBJECTS_ON_QUAD_CELL(K, I, J);

//                COEFF_SEC += WEIGHT * FUN * PHI(J);
//            }

//            m_duos.push_back ({ID(J), COEFF_SEC});
//        }

//    }

//    if (view) {ENDFUN;}

//    return;
//}

///** ********************************************************************** *
//  * ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY
//  * ********************************************************************** *
//  *
//  * U_D = U on boundary
//  * [C] = coeff penalization
//  *
//  * MATRIX :
//  * -1 <PHI(J), GRAD_PHI(I).NORMAL>
//  * -1 <GRAD_PHI(J).NORMAL | PHI(I)>
//  * +1 <[C] * PHI(J) | PHI(I)>
//  *
//  * SECOND MEMBRE :
//  * -1 <GRAD_PHI(J).NORMAL | U_D>
//  * +1 <[C] * PHI(J) | U_D>
//  **/

//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();
//    int numEdges  = (*m_mesh)->GetNumberOfEdges ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    if (*m_collect)
//    {
//        m_applicableCells.clear ();
//        m_applicableEdges.clear ();

//        m_applicableCells.resize (static_cast<std::size_t>(numCells), 0);
//        m_applicableEdges.resize (static_cast<std::size_t>(numEdges), 0);
//    }

//    double penal = *m_coeff_pen / ((*m_mesh)->GetPrescribedSize () * 1.);

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}

//        if (*m_collect)
//            m_applicableCells.at (static_cast<std::size_t>(CELLID)) += 1;

//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                for (std::size_t EDGEID = 0; EDGEID < numOfEdgesOnCell; ++EDGEID)
//                {
//                    OBJECTS_ON_EDGE(EDGEID);

//                    if (tagPhysicalEdge != m_tagToApply)
//                        continue;

//                    if (*m_collect)
//                        m_applicableEdges.at (static_cast<std::size_t>(edge->GetGlobalIndex ())) += 1;

//                    double COEFF_MATRIX = 0.;
//                    double COEFF_SEC    = 0.;

//                    for (std::size_t K = 0; K < QUAD4EDGE->npts; ++K)
//                    {
//                        OBJECTS_ON_QUAD_EDGE(K, I, J);

//                        //!! -1 <PHI(J), GRAD_PHI(I).NORMAL>
//                        COEFF_MATRIX -= WEIGHT * PHI(J) * (GRAD_PHI(I) | NORMAL);
//                        //!! -1 <GRAD_PHI(J).NORMAL | PHI(I)>
//                        COEFF_MATRIX -= WEIGHT * (GRAD_PHI(J) | NORMAL) * PHI(I);
//                        //!! +1 <[C] * PHI(J) | PHI(I)>
//                        COEFF_MATRIX += WEIGHT * penal * PHI(J) * PHI(I);

//                        if (I == 0)
//                        {
//                            //!! -1 <GRAD_PHI(J).NORMAL | U_D>
//                            COEFF_SEC -= WEIGHT * (GRAD_PHI(J) | NORMAL) * FUN;
//                            //!! +1 <[C] * PHI(J) | U_D>
//                            COEFF_SEC += WEIGHT * penal * PHI(J) * FUN;
//                        }

//                    }

//                    m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});
//                    if (I == 0)
//                        m_duos.push_back ({ID(J), COEFF_SEC});
//                }
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    //    ERROR << "ITEM DIRICHLET " << m_tagToApply << ENDLINE;
//    //    std::cout << Eigen::MatrixXd(*this->CastSparseMatrix ()) << ENDLINE;
//    return;
//}

///** ********************************************************************** *
// * ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY
// * *********************************************************************** *
// *
// * t_n = du/dn (neumann condition)
// *
// * MATRIX :
// * +1 <PHI(J), GRAD_PHI(I).NORMAL>
// *
// * SECOND MEMBRE :
// * + <PHI(J), t_n>
// **/

//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}


//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                for (std::size_t EDGEID = 0; EDGEID < numOfEdgesOnCell; ++EDGEID)
//                {
//                    OBJECTS_ON_EDGE (EDGEID);

//                    if (tagPhysicalEdge != m_tagToApply)
//                        continue;

//                    double COEFF_MATRIX = 0.;
//                    double COEFF_SEC    = 0.;

//                    for (std::size_t K = 0; K < QUAD4EDGE->npts; ++K)
//                    {
//                        OBJECTS_ON_QUAD_EDGE(K, I, J);

//                        COEFF_MATRIX += WEIGHT * PHI(J) * (GRAD_PHI(I) | NORMAL);
//                        if (I == 0)
//                            COEFF_SEC += WEIGHT * PHI(J) * FUN;
//                    }

//                    m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});
//                    if (I == 0)
//                        m_duos.push_back ({ID(J), COEFF_SEC});
//                }
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    return;
//}

///** ********************************************************************** *
// * ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM
// * *********************************************************************** *
// *
// * DISP_X = PHI_X + GRAD_PHI_X.d = displacement vector
// * [C] = coeff penalization
// * U_D = U on boundary
// * {X} = surrogate resp of X
// *
// * MATRIX :
// * -1 <PHI(J) | GRAD_PHI(I).{NORMAL}>
// * -1 <GRAD_PHI(J).{NORMAL} | DISP_I>
// * +1 <[C] * DISP_J | DISP_I>
// *
// * SECOND MEMBRE :
// * -1 <GRAD_PHI(J).{NORMAL} | {U_D}>
// * +1 <[C] * DEP_J | {U_D}>
// **/

//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}


//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                for (std::size_t EDGEID = 0; EDGEID < numOfEdgesOnCell; ++EDGEID)
//                {
//                    OBJECTS_ON_EDGE (EDGEID);

//                    if (tagPhysicalEdge != m_tagToApply)
//                        continue;

//                    double COEFF_MATRIX = 0.;
//                    double COEFF_SEC    = 0.;

//                    for (std::size_t K = 0; K < QUAD4EDGE->npts; ++K)
//                    {
//                        OBJECTS_ON_QUAD_EDGE(K, I, J);

//                        Point d;
//                        GetDisplacementVectorAtPoint (&POINT_REAL, m_emb, &d);
//                        d = locCell.Jac_cmp * locEdge.Jac_cmp * d;

//                        COEFF_MATRIX -= WEIGHT * PHI(J) * (GRAD_PHI(I) | NORMAL);
//                        COEFF_MATRIX -= WEIGHT * (GRAD_PHI(J) | NORMAL) * (PHI(I) + (GRAD_PHI(I) | d));
//                        COEFF_MATRIX += WEIGHT * *m_coeff_pen * (PHI(J) + (GRAD_PHI(J) | d)) * (PHI(I) + (GRAD_PHI(I) | d));

//                        if (I == 0)
//                        {
//                            COEFF_SEC -= WEIGHT * (GRAD_PHI(J) | NORMAL) * FUN;
//                            COEFF_SEC += WEIGHT * (PHI(J) + (GRAD_PHI(J) | d)) * FUN;
//                        }
//                    }

//                    m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});
//                    if (I == 0)
//                        m_duos.push_back ({ID(J), COEFF_SEC});
//                }
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    return;
//}

///** ********************************************************************** *
// * ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY_SBM
// * ********************************************************************** **/

//template<>
//void ItemSolver<ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY_SBM>::Compute (double time, bool view, bool force)
//{
//    VOID_USE(time); VOID_USE(view); VOID_USE(force);

//    if (!(m_variableOverTime || force))
//        return;

//    m_duos.clear ();
//    m_triplets.clear ();

//    int numCells  = (*m_mesh)->GetNumberOfCells ();

//    m_duos.reserve (static_cast<std::size_t>(32 * numCells));
//    m_triplets.reserve (static_cast<std::size_t>(32 * numCells));

//    for (int CELLID = 0; CELLID < numCells; ++CELLID)
//    {
//        OBJECTS_ON_CELL(CELLID);

//        if (view) {DISPLAY_PERCENTS(CELLID, numCells);}


//        for (std::size_t I = 0; I < numOfPointsOnCell; ++I)
//        {
//            OBJECTS_ON_POINT(I);

//            for (std::size_t J = 0; J < numOfPointsOnCell; ++J)
//            {
//                OBJECTS_ON_POINT(J);

//                double COEFF_MATRIX = 0.;
//                double COEFF_SEC    = 0.;

//                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)

//                {
//                    OBJECTS_ON_QUAD_CELL(K, I, J);

//                    COEFF_MATRIX += WEIGHT * PHI(J) * PHI(I);

//                    if (I == 0)
//                        COEFF_SEC += WEIGHT * 0.;
//                }

//                m_triplets.push_back ({ID(I), ID(J), COEFF_MATRIX});

//                if (I == 0)
//                    m_duos.push_back ({ID(J), COEFF_SEC});
//            }
//        }
//    }

//    if (view)
//        ENDFUN;

//    return;
//}

//#include "undef_macros_item.h"
