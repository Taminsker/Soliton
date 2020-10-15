#include <Solver/SuperSolver/SuperItem/superitem.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>
#include <Algorithms/algorithms.h>

#include "def_macros_items.h"

template<>
void
SuperItem::InternalCompute<ITEM_T::EMPTY> ()
{
    SuperItem *ITEM = this;
    VOID_USE(ITEM);

    return;
}

/** ********************************************************************** *
  * ITEM_SOLVER_TYPE::PHI_PHI
  * ********************************************************************** *
  *
  * +<PHI(J), PHI(I)>
  *
  **/
template<>
void
SuperItem::InternalCompute<ITEM_T::PHI_PHI> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}

        if (TAG_PHYSICAL_CELL != ITEM->GetTag2Apply ())
            continue;

        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                double COEFF_MATRIX = 0.;
                double COEFF_SEC    = 0.;

                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
                {
                    OBJECTS_ON_QUAD_CELL(K, I, J);

                    COEFF_MATRIX += WEIGHT * PHI(J) * PHI(I);
                    if (I == 0)
                        COEFF_SEC += WEIGHT * 0.;
                }

                TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});

                if (I == 0)
                    DUOS.push_back ({ID(J), COEFF_SEC});
            }
        }
    }

    ITEM_FINALIZE;
}

/** ********************************************************************** *
  * ITEM_SOLVER_TYPE::DAMPING_PHI_PHI
  * ********************************************************************** *
  *  +<PHI(J), chi PHI(I)>
  **/
template<>
void
SuperItem::InternalCompute<ITEM_T::DAMPING_PHI_PHI> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

    auto FunDamp = [=](Point* p) -> double
    {return GetDampingCoeffFor (p, MESH, ITEM->GetTag2Apply (), ITEM->GetCoeffPen ());};


    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}


        if (TAG_PHYSICAL_CELL != ITEM->GetTag2Apply ())
            continue;

        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                double COEFF_MATRIX = 0.;

                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
                {
                    OBJECTS_ON_QUAD_CELL(K, I, J);

                    COEFF_MATRIX += WEIGHT *  1. / FunDamp (&POINT_REAL) * PHI(J) * PHI(I);
                }

                TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});
            }
        }
    }

    ITEM_FINALIZE;
}

/** ********************************************************************** *
  * ITEM_SOLVER_TYPE::GRAD_GRAD
  * ********************************************************************** *
  *  +<GRAD_PHI(J), GRAD_PHI(I)>
  **/
template<>
void
SuperItem::InternalCompute<ITEM_T::GRAD_GRAD> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

//    if (*m_collect)
//    {
//        m_applicableCells.clear ();
//        m_applicableEdges.clear ();

//        m_applicableCells.resize (static_cast<std::size_t>(numCells), 0x0);
//        m_applicableEdges.resize (static_cast<std::size_t>(numEdges), 0x0);
//    }

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}

        if (TAG_PHYSICAL_CELL != ITEM->GetTag2Apply ())
            continue;

//        if (*m_collect)
//            m_applicableCells.at (static_cast<std::size_t>(CELLID)) += 1;

        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                double COEFF_MATRIX = 0.;

                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
                {
                    OBJECTS_ON_QUAD_CELL(K, I, J);

                    COEFF_MATRIX += WEIGHT * (GRAD_PHI(J) | GRAD_PHI(I));
                }

                TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});
            }
        }
    }

    ITEM_FINALIZE;
}

/** ********************************************************************** *
  * ITEM_SOLVER_TYPE::SECOND_MEMBER
  * ********************************************************************** *
  *  +<F, PHI(J)>
  **/
template<>
void
SuperItem::InternalCompute<ITEM_T::SECOND_MEMBER> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

//    if (*m_collect)
//    {
//        m_applicableCells.resize (static_cast<std::size_t>(numCells), 0);
//        m_applicableEdges.resize (static_cast<std::size_t>(numEdges), 0);
//    }

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}

        if (TAG_PHYSICAL_CELL != ITEM->GetTag2Apply ())
            continue;

//        if (*m_collect)
//            m_applicableCells.at (static_cast<std::size_t>(CELLID)) += 1;

        std::size_t I = 0;
        OBJECTS_ON_POINT(I);

        for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
        {
            OBJECTS_ON_POINT(J);

            double COEFF_SEC    = 0.;

            for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)
            {
                OBJECTS_ON_QUAD_CELL(K, I, J);

                COEFF_SEC += WEIGHT * FUN * PHI(J);
            }

            DUOS.push_back ({ID(J), COEFF_SEC});
        }

    }

    ITEM_FINALIZE;
}

/** ********************************************************************** *
  * ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY
  * ********************************************************************** *
  *
  * U_D = U on boundary
  * [C] = coeff penalization
  *
  * MATRIX :
  * -1 <PHI(J), GRAD_PHI(I).NORMAL>
  * -1 <GRAD_PHI(J).NORMAL | PHI(I)>
  * +1 <[C] * PHI(J) | PHI(I)>
  *
  * SECOND MEMBRE :
  * -1 <GRAD_PHI(J).NORMAL | U_D>
  * +1 <[C] * PHI(J) | U_D>
  **/

template<>
void
SuperItem::InternalCompute<ITEM_T::DIRICHLET_BOUNDARY> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

//    if (*m_collect)
//    {
//        m_applicableCells.clear ();
//        m_applicableEdges.clear ();

//        m_applicableCells.resize (static_cast<std::size_t>(numCells), 0);
//        m_applicableEdges.resize (static_cast<std::size_t>(numEdges), 0);
//    }

    double penal = ITEM->GetCoeffPen () / (MESH->GetPrescribedSize () * 1.);

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}

//        if (*m_collect)
//            m_applicableCells.at (static_cast<std::size_t>(CELLID)) += 1;

        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                for (std::size_t EDGEID = 0; EDGEID < NUM_EDGES_ON_CELL; ++EDGEID)
                {
                    OBJECTS_ON_EDGE(EDGEID);

                    if (TAG_PHYSICAL_EDGE != ITEM->GetTag2Apply ())
                        continue;

//                    if (*m_collect)
//                        m_applicableEdges.at (static_cast<std::size_t>(edge->GetGlobalIndex ())) += 1;

                    double COEFF_MATRIX = 0.;
                    double COEFF_SEC    = 0.;

                    for (std::size_t K = 0; K < QUAD4EDGE->npts; ++K)
                    {
                        OBJECTS_ON_QUAD_EDGE(K, I, J);

                        //!! -1 <PHI(J), GRAD_PHI(I).NORMAL>
                        COEFF_MATRIX -= WEIGHT * PHI(J) * (GRAD_PHI(I) | NORMAL);
                        //!! -1 <GRAD_PHI(J).NORMAL | PHI(I)>
                        COEFF_MATRIX -= WEIGHT * (GRAD_PHI(J) | NORMAL) * PHI(I);
                        //!! +1 <[C] * PHI(J) | PHI(I)>
                        COEFF_MATRIX += WEIGHT * penal * PHI(J) * PHI(I);

                        if (I == 0)
                        {
                            //!! -1 <GRAD_PHI(J).NORMAL | U_D>
                            COEFF_SEC -= WEIGHT * (GRAD_PHI(J) | NORMAL) * FUN;
                            //!! +1 <[C] * PHI(J) | U_D>
                            COEFF_SEC += WEIGHT * penal * PHI(J) * FUN;
                        }

                    }

                    TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});
                    if (I == 0)
                        DUOS.push_back ({ID(J), COEFF_SEC});
                }
            }
        }
    }

   ITEM_FINALIZE;
}

/** ********************************************************************** *
 * ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY
 * *********************************************************************** *
 *
 * t_n = du/dn (neumann condition)
 *
 * MATRIX :
 * +1 <PHI(J), GRAD_PHI(I).NORMAL>
 *
 * SECOND MEMBRE :
 * + <PHI(J), t_n>
 **/

template<>
void
SuperItem::InternalCompute<ITEM_T::NEUMANN_BOUNDARY> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}

        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                for (std::size_t EDGEID = 0; EDGEID < NUM_EDGES_ON_CELL; ++EDGEID)
                {
                    OBJECTS_ON_EDGE (EDGEID);

                    if (TAG_PHYSICAL_EDGE != ITEM->GetTag2Apply ())
                        continue;

                    double COEFF_MATRIX = 0.;
                    double COEFF_SEC    = 0.;

                    for (std::size_t K = 0; K < QUAD4EDGE->npts; ++K)
                    {
                        OBJECTS_ON_QUAD_EDGE(K, I, J);

                        COEFF_MATRIX += WEIGHT * PHI(J) * (GRAD_PHI(I) | NORMAL);
                        if (I == 0)
                            COEFF_SEC += WEIGHT * PHI(J) * FUN;
                    }

                    TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});
                    if (I == 0)
                        DUOS.push_back ({ID(J), COEFF_SEC});
                }
            }
        }
    }

    ITEM_FINALIZE;
}

/** ********************************************************************** *
 * ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM
 * *********************************************************************** *
 *
 * DISP_X = PHI_X + GRAD_PHI_X.d = displacement vector
 * [C] = coeff penalization
 * U_D = U on boundary
 * {X} = surrogate resp of X
 *
 * MATRIX :
 * -1 <PHI(J) | GRAD_PHI(I).{NORMAL}>
 * -1 <GRAD_PHI(J).{NORMAL} | DISP_I>
 * +1 <[C] * DISP_J | DISP_I>
 *
 * SECOND MEMBRE :
 * -1 <GRAD_PHI(J).{NORMAL} | {U_D}>
 * +1 <[C] * DEP_J | {U_D}>
 **/

template<>
void
SuperItem::InternalCompute<ITEM_T::DIRICHLET_BOUNDARY_SBM> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

    double penal = ITEM->GetCoeffPen () / (MESH->GetPrescribedSize () * 1.);

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}

        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                for (std::size_t EDGEID = 0; EDGEID < NUM_EDGES_ON_CELL; ++EDGEID)
                {
                    OBJECTS_ON_EDGE (EDGEID);

                    if (TAG_PHYSICAL_EDGE != ITEM->GetTag2Apply ())
                        continue;

                    double COEFF_MATRIX = 0.;
                    double COEFF_SEC    = 0.;

                    for (std::size_t K = 0; K < QUAD4EDGE->npts; ++K)
                    {
                        OBJECTS_ON_QUAD_EDGE(K, I, J);

                        Point d;
                        GetDisplacementVectorAtPoint (&POINT_REAL, ITEM->GetTargetMesh (), &d);

                        d = locCell.Jac_cmp * locEdge.Jac_cmp * d;

                        COEFF_MATRIX -= WEIGHT * PHI(J) * (GRAD_PHI(I) | NORMAL);
                        COEFF_MATRIX -= WEIGHT * (GRAD_PHI(J) | NORMAL) * (PHI(I) + (GRAD_PHI(I) | d));
                        COEFF_MATRIX += WEIGHT * penal * (PHI(J) + (GRAD_PHI(J) | d)) * (PHI(I) + (GRAD_PHI(I) | d));

                        if (I == 0)
                        {
                            COEFF_SEC -= WEIGHT * (GRAD_PHI(J) | NORMAL) * FUN;
                            COEFF_SEC += WEIGHT * (PHI(J) + (GRAD_PHI(J) | d)) * FUN;
                        }
                    }

                    TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});
                    if (I == 0)
                        DUOS.push_back ({ID(J), COEFF_SEC});
                }
            }
        }
    }

    ITEM_FINALIZE;
}

/** ********************************************************************** *
 * ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY_SBM
 * ********************************************************************** **/

template<>
void
SuperItem::InternalCompute<ITEM_T::NEUMANN_BOUNDARY_SBM> ()
{
    SuperItem *ITEM = this;

    ITEM_BEGIN;

    for (int CELLID = 0; CELLID < NUMCELLS; ++CELLID)
    {
        OBJECTS_ON_CELL(CELLID);

        if (ITEM->GetView ()) {DISPLAY_PERCENTS(CELLID, NUMCELLS);}


        for (std::size_t I = 0; I < NUM_POINTS_ON_CELL; ++I)
        {
            OBJECTS_ON_POINT(I);

            for (std::size_t J = 0; J < NUM_POINTS_ON_CELL; ++J)
            {
                OBJECTS_ON_POINT(J);

                double COEFF_MATRIX = 0.;
                double COEFF_SEC    = 0.;

                for (std::size_t K = 0; K < QUAD4CELL->npts; ++K)

                {
                    OBJECTS_ON_QUAD_CELL(K, I, J);

                    COEFF_MATRIX += WEIGHT * PHI(J) * PHI(I);

                    if (I == 0)
                        COEFF_SEC += WEIGHT * 0.;
                }

                TRIPLETS.push_back ({ID(I), ID(J), COEFF_MATRIX});

                if (I == 0)
                    DUOS.push_back ({ID(J), COEFF_SEC});
            }
        }
    }

    ITEM_FINALIZE;
}

#include "undef_macros_item.h"
