#ifndef DEF_MACROS_ITEMS_H
#define DEF_MACROS_ITEMS_H

#include <ProgDef/proddef.h>

#define CELLID              MAKE_SOLITON_RESERVE(CellId)
#define EDGEID              MAKE_SOLITON_RESERVE(EdgeId)
#define I                   MAKE_SOLITON_RESERVE(i)
#define J                   MAKE_SOLITON_RESERVE(j)
#define K                   MAKE_SOLITON_RESERVE(k)

#define ID(X)               MAKE_SOLITON_RESERVE(CONCAT_STR_2(id_pt_, X))
#define FUN_PHI(X)          MAKE_SOLITON_RESERVE(CONCAT_STR_2(fun_phi_, X))
#define FUN_GRAD_PHI(X)     MAKE_SOLITON_RESERVE(CONCAT_STR_2(fun_grad_phi_, X))
#define PHI(X)              MAKE_SOLITON_RESERVE(CONCAT_STR_2(phi_, X))
#define GRAD_PHI(X)         MAKE_SOLITON_RESERVE(CONCAT_STR_2(grad_phi_, X))

#define FUN                 MAKE_SOLITON_RESERVE(val_f)
#define NORMAL              MAKE_SOLITON_RESERVE(normal)
#define WEIGHT              MAKE_SOLITON_RESERVE(w)
#define POINT_TRAN          MAKE_SOLITON_RESERVE(pt_tran)
#define POINT_REAL          MAKE_SOLITON_RESERVE(pt_real)
#define POINT_INT           MAKE_SOLITON_RESERVE(pt)

#define QUAD4CELL           MAKE_SOLITON_RESERVE(quad4Cell)
#define QUAD4EDGE           MAKE_SOLITON_RESERVE(quad4Edge)

#define COEFF_MATRIX        MAKE_SOLITON_RESERVE(coeff_matrix)
#define COEFF_SEC           MAKE_SOLITON_RESERVE(coeff_second_member)


#define DISPLAY_PERCENTS(X) \
        int percent = static_cast<int>(static_cast<double>(X + 1) / static_cast<double>((*m_mesh)->GetNumberOfCells ()) * 100.); \
        std::cout << "\r" << COLOR_BLUE << "[" << percent << "%]  " << COLOR_DEFAULT << *this << std::flush

#define OBJECTS_ON_CELL(X) \
        Cell* cell = (*m_mesh)->GetCell (X); \
        FEBase *febase4Cell = m_festore->GetElementFor (cell); \
        std::vector<Edge*> *listEdgesOnCell = cell->GetEdges (); \
        std::size_t numOfEdgesOnCell = listEdgesOnCell->size (); \
        QuadStore::QuadObject *QUAD4CELL = m_quadstore->Get (febase4Cell); \
        TAG_PHYSICAL tagPhysicalCell = TAG_PHYSICAL((*m_mesh)->GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL)->vec.at (std::size_t(cell->GetGlobalIndex ()))); \
        std::size_t numOfPointsOnCell = std::size_t(cell->GetNumberOfPoints()); \
        (void) QUAD4CELL; \
        (void) numOfEdgesOnCell; \
        (void) tagPhysicalCell; \
        (void) numOfPointsOnCell

#define OBJECTS_ON_POINT(X) \
        LambdaOnPoint<double> FUN_PHI(X) = febase4Cell->GetPhi (X); \
        LambdaOnPoint<Point> FUN_GRAD_PHI(X) = febase4Cell->GetGradPhi (X); \
        int ID(X) = cell->GetPoints ()->at (X)->GetGlobalIndex (); \
        (void) ID(X)

#define OBJECTS_ON_EDGE(X) \
        Edge*                   edge            = cell->GetEdges()->at(X); \
        Edge*                   edgeOnCellRef   = febase4Cell->GetEdge (X); \
        FEBase*                 febase4Edge     = m_festore->GetElementFor (edgeOnCellRef); \
        QuadStore::QuadObject*  QUAD4EDGE       = m_quadstore->Get (febase4Edge); \
        TAG_PHYSICAL            tagPhysicalEdge  = TAG_PHYSICAL((*m_mesh)->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL)->vec.at (std::size_t(edge->GetGlobalIndex ()))); \
        (void) tagPhysicalEdge

#define OBJECTS_ON_QUAD_CELL(V_K, V_I, V_J) \
        Point   POINT_INT = QUAD4CELL->pts [V_K]; \
        FELocalInfos locCell; febase4Cell->LocalCompute (&POINT_INT, &locCell); \
        Point POINT_REAL = febase4Cell->TransformRefToEle (&POINT_INT); \
        double WEIGHT = locCell.detJac * QUAD4CELL->w [V_K]; \
        Point GRAD_PHI(V_I) = locCell.JacInvT * FUN_GRAD_PHI(V_I) (&POINT_INT); \
        Point GRAD_PHI(V_J) = locCell.JacInvT * FUN_GRAD_PHI(V_J) (&POINT_INT); \
        double PHI(V_I) = FUN_PHI(V_I) (&POINT_INT); \
        double PHI(V_J) = FUN_PHI(V_J) (&POINT_INT); \
        double FUN = m_fun(POINT_REAL, time); \
        (void) GRAD_PHI(V_I); (void) GRAD_PHI(V_J); (void) PHI(V_I); (void) PHI(V_J); (void) FUN



#define OBJECTS_ON_QUAD_EDGE(V_K, V_I, V_J) \
        Point POINT_INT = QUAD4EDGE->pts [V_K]; \
        Point POINT_TRAN = febase4Edge->TransformRefToEle (&POINT_INT); \
        Point POINT_REAL = febase4Cell->TransformRefToEle (&POINT_TRAN); \
        FELocalInfos locEdge;  febase4Edge->LocalCompute (&POINT_INT, &locEdge); \
        FELocalInfos locCell; febase4Cell->LocalCompute (&POINT_TRAN, &locCell); \
        double WEIGHT = locEdge.detJac * locCell.detJac * QUAD4EDGE->w [V_K]; \
        Point GRAD_PHI(V_I) = locCell.JacInvT * locEdge.JacInvT * FUN_GRAD_PHI(V_I) (&POINT_TRAN); \
        Point GRAD_PHI(V_J) = locCell.JacInvT * locEdge.JacInvT * FUN_GRAD_PHI(V_J) (&POINT_TRAN); \
        double PHI(V_I) = FUN_PHI(V_I) (&POINT_TRAN);\
        double PHI(V_J) = FUN_PHI(J) (&POINT_TRAN);\
        double FUN = m_fun(POINT_REAL, time);\
        Point NORMAL = locEdge.normalCell_ref; \
        (void) GRAD_PHI(V_I); (void) GRAD_PHI(V_J); (void) PHI(V_I); (void) PHI(V_J); (void) FUN; (void) NORMAL

#endif // DEF_MACROS_ITEMS_H
