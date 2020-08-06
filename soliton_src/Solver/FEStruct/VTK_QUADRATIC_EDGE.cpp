#include "festruct.h"


//template <>
//void LagrangeFE<VTK_QUADRATIC_EDGE>::Build ()
//{
////    npts = 3;
////    coor = {{-1}, {0}, {1}};
////    volRef = 2;

////    fun.resize (npts);
////    der.resize (npts);

////    fun [0]   = LAMBDA_EVAL (-2 * p.x * (1 - p.x));
////    fun [1]   = LAMBDA_EVAL ( 2 * (1 - p.x * p.x));
////    fun [2]   = LAMBDA_EVAL ( 2 * p.x * (1 + p.x));

////    der [0].x = LAMBDA_EVAL (2 * (-1 + 2 * p.x));
////    der [1].x = LAMBDA_EVAL (2 * (-4 *  p.x));
////    der [2].x = LAMBDA_EVAL (2 * (1 + 2 * p.x));

//    return;
//}

//template <>
//Matrix LagrangeFE<VTK_QUADRATIC_EDGE>::Jac (std::vector <Point*> P)
//{
//    if (P.size() != 3)
//        return {{0}};

//    Matrix mat = {{0.5 * (P[1]->x - P[0]->x)}};

//    return mat;
//}
