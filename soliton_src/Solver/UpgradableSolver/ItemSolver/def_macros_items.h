//#ifndef DEF_MACROS_ITEMS_H
//#define DEF_MACROS_ITEMS_H

//#include <ProgDef/proddef.h>

//#define CELLID              MAKE_SOLITON_RESERVE(CellId)
//#define EDGEID              MAKE_SOLITON_RESERVE(EdgeId)
//#define I                   MAKE_SOLITON_RESERVE(i)
//#define J                   MAKE_SOLITON_RESERVE(j)
//#define K                   MAKE_SOLITON_RESERVE(k)

//#define ID(X)               MAKE_SOLITON_RESERVE(CONCAT_STR_2(id_pt_, X))
//#define FUN_PHI(X)          MAKE_SOLITON_RESERVE(CONCAT_STR_2(fun_phi_, X))
//#define FUN_GRAD_PHI(X)     MAKE_SOLITON_RESERVE(CONCAT_STR_2(fun_grad_phi_, X))
//#define PHI(X)              MAKE_SOLITON_RESERVE(CONCAT_STR_2(phi_, X))
//#define GRAD_PHI(X)         MAKE_SOLITON_RESERVE(CONCAT_STR_2(grad_phi_, X))

//#define FUN                 MAKE_SOLITON_RESERVE(val_f)
//#define NORMAL              MAKE_SOLITON_RESERVE(normal)
//#define WEIGHT              MAKE_SOLITON_RESERVE(w)
//#define POINT_TRAN          MAKE_SOLITON_RESERVE(pt_tran)
//#define POINT_REAL          MAKE_SOLITON_RESERVE(pt_real)
//#define POINT_INT           MAKE_SOLITON_RESERVE(pt)

//#define QUAD4CELL           MAKE_SOLITON_RESERVE(quad4Cell)
//#define QUAD4EDGE           MAKE_SOLITON_RESERVE(quad4Edge)

//#define COEFF_MATRIX        MAKE_SOLITON_RESERVE(coeff_matrix)
//#define COEFF_SEC           MAKE_SOLITON_RESERVE(coeff_second_member)


//#define DISPLAY_PERCENTS(X, Y) \
//        int percent = static_cast<int>(static_cast<double>(X + 1) / static_cast<double>(Y) * 100.); \
//        std::cout << "\r\t";  \
//        TREE_BRANCH << COLOR_BLUE << "[" << percent << "%] " << COLOR_DEFAULT << *this << " " << COLOR_MAGENTA << "[Cell " << COLOR_BLUE << X+1 << "/" << Y << COLOR_MAGENTA << " " << to_string(cell->GetTypeVTK()) << "]                         " << std::flush

//#define OBJECTS_ON_CELL(X) \
//        Cell* cell = (*m_mesh)->GetCell (X); \
//        FEBase *febase4Cell = m_festore->GetElementFor (cell); \
//        std::size_t numOfEdgesOnCell = cell->GetEdges ()->size (); \
//        QuadStore::QuadObject *QUAD4CELL = m_quadstore->Get (febase4Cell); \
//        PHYS tagPhysicalCell = PHYS((*m_mesh)->GetCellsData ()->GetIntArrays ()->Get (NAME_PHYS)->vec.at (std::size_t(cell->GetGlobalIndex ()))); \
//        std::size_t numOfPointsOnCell = std::size_t(cell->GetNumberOfPoints()); \
//    \
//        VOID_USE(QUAD4CELL); \
//        VOID_USE(numOfEdgesOnCell); \
//        VOID_USE(tagPhysicalCell); \
//        VOID_USE(numOfPointsOnCell)

//#define OBJECTS_ON_POINT(X) \
//        LambdaOnPoint<double> FUN_PHI(X) = febase4Cell->GetPhi (X); \
//        LambdaOnPoint<Point> FUN_GRAD_PHI(X) = febase4Cell->GetGradPhi (X); \
//        int ID(X) = cell->GetPoints ()->at (X)->GetGlobalIndex (); \
//        VOID_USE(ID(X))

//#define OBJECTS_ON_EDGE(X) \
//        Edge*                   edge            = cell->GetEdges()->at(X); \
//        Edge*                   edgeOnCellRef   = febase4Cell->GetEdge (X); \
//        FEBase*                 febase4Edge     = m_festore->GetElementFor (edgeOnCellRef); \
//        QuadStore::QuadObject*  QUAD4EDGE       = m_quadstore->Get (febase4Edge); \
//        PHYS            tagPhysicalEdge  = PHYS((*m_mesh)->GetEdgesData ()->GetIntArrays ()->Get (NAME_PHYS)->vec.at (std::size_t(edge->GetGlobalIndex ()))); \
//        VOID_USE(tagPhysicalEdge)

//#define OBJECTS_ON_QUAD_CELL(V_K, V_I, V_J) \
//        Point   POINT_INT = QUAD4CELL->pts [V_K]; \
//        Point POINT_REAL = febase4Cell->TransformRefToEle (&POINT_INT); \
//        FELocalInfos locCell; febase4Cell->LocalCompute (&POINT_INT, &locCell); \
//    \
//        double WEIGHT = locCell.detJac * QUAD4CELL->w [V_K]; \
//    \
//        Point GRAD_PHI(V_I) = locCell.JacInvT * FUN_GRAD_PHI(V_I) (&POINT_INT); \
//        Point GRAD_PHI(V_J) = locCell.JacInvT * FUN_GRAD_PHI(V_J) (&POINT_INT); \
//    \
//        double PHI(V_I) = FUN_PHI(V_I) (&POINT_INT); \
//        double PHI(V_J) = FUN_PHI(V_J) (&POINT_INT); \
//    \
//        double FUN = m_fun(POINT_REAL, time); \
//    \
//        VOID_USE(GRAD_PHI(V_I)); VOID_USE(GRAD_PHI(V_J)); \
//        VOID_USE(PHI(V_I)); VOID_USE(PHI(V_J)); \
//        VOID_USE(FUN); VOID_USE(WEIGHT)

//#define OBJECTS_ON_QUAD_EDGE(V_K, V_I, V_J) \
//        Point POINT_INT = QUAD4EDGE->pts [V_K]; \
//        Point POINT_TRAN = febase4Edge->TransformRefToEle (&POINT_INT); \
//        Point POINT_REAL = febase4Cell->TransformRefToEle (&POINT_TRAN); \
//    \
//        FELocalInfos locEdge; febase4Edge->LocalCompute (&POINT_INT, &locEdge); \
//        FELocalInfos locCell; febase4Cell->LocalCompute (&POINT_TRAN, &locCell); \
//    \
//        double WEIGHT = locEdge.detJac * locCell.detJac * QUAD4EDGE->w [V_K]; \
//    \
//        Point GRAD_PHI(V_I) = locEdge.JacInvT * locCell.JacInvT * FUN_GRAD_PHI(V_I) (&POINT_TRAN); \
//        Point GRAD_PHI(V_J) = locEdge.JacInvT * locCell.JacInvT * FUN_GRAD_PHI(V_J) (&POINT_TRAN); \
//    \
//        double PHI(V_I) = FUN_PHI(V_I) (&POINT_TRAN);\
//        double PHI(V_J) = FUN_PHI(V_J) (&POINT_TRAN);\
//    \
//        double FUN = m_fun(POINT_REAL, time); \
//        Point NORMAL = locEdge.normalCell_ref; \
//    \
//        VOID_USE(GRAD_PHI(V_I)); VOID_USE(GRAD_PHI(V_J)); \
//        VOID_USE(PHI(V_I)); VOID_USE(PHI(V_J)); \
//        VOID_USE(FUN); VOID_USE(NORMAL); VOID_USE(WEIGHT)


///*
//#define OBJECTS_ON_QUAD_EDGE(V_K, V_I, V_J) \
//        Point POINT_INT = QUAD4EDGE->pts [V_K]; \
//        Point POINT_TRAN = febase4Edge->TransformRefToEle (&POINT_INT); \
//        Point POINT_REAL = febase4Cell->TransformRefToEle (&POINT_TRAN); \
//    \
//        FELocalInfos locEdge; febase4Edge->LocalCompute (&POINT_INT, &locEdge); \
//        FELocalInfos locCell; febase4Cell->LocalCompute (&POINT_TRAN, &locCell); \
//    \
//        double WEIGHT = locCell.detJac * QUAD4EDGE->w [V_K]; \
//    \
//        Point GRAD_PHI(V_I) = locCell.JacInvT_cmp * FUN_GRAD_PHI(V_I) (&POINT_TRAN); \
//        Point GRAD_PHI(V_J) = locCell.JacInvT_cmp * FUN_GRAD_PHI(V_J) (&POINT_TRAN); \
//    \
//        double PHI(V_I) = FUN_PHI(V_I) (&POINT_TRAN);\
//        double PHI(V_J) = FUN_PHI(V_J) (&POINT_TRAN);\
//    \
//        double FUN = m_fun(POINT_REAL, time);\
//        Point NORMAL = locEdge.normalEdge_ref; \
//    \
//        VOID_USE(GRAD_PHI(V_I)); VOID_USE(GRAD_PHI(V_J)); \
//        VOID_USE(PHI(V_I)); VOID_USE(PHI(V_J)); \
//        VOID_USE(FUN); VOID_USE(NORMAL)
//*/

//#endif // DEF_MACROS_ITEMS_H
