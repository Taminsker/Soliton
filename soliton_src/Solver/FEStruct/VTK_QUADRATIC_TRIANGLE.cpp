#include "festruct.h"

template <>
void LagrangeFE<VTK_QUADRATIC_TRIANGLE>::Build ()
{
//    npts = 6;
//    coor = {{0, 0}, {0.5, 0}, {1, 0}, {0.5, 0.5}, {0, 1}, {0, 0.5}};

//    fun.resize (npts);
//    der.resize (npts);

//    fun [0] = LAMBDA_EVAL(-(1. - p.x - p.y) * (1 - 2 * (1. - p.x - p.y)););
//    fun [1] = LAMBDA_EVAL(4 * p.x * (1. - p.x - p.y););
//    fun [2] = LAMBDA_EVAL(- p.x * (1 - 2 * p.x););
//    fun [3] = LAMBDA_EVAL(4 * p.x * p.y;);
//    fun [4] = LAMBDA_EVAL(- p.y * (1 - 2 * p.y););
//    fun [5] = LAMBDA_EVAL(4 * p.y * (1. - p.x - p.y););

//    der [0].x = LAMBDA_EVAL (1 - 4 * (1. - p.x - p.y));
//    der [0].y = LAMBDA_EVAL (1 - 4 * (1. - p.x - p.y));

//    der [1].x = LAMBDA_EVAL (4 * ((1. - p.x - p.y) - p.x));
//    der [1].y = LAMBDA_EVAL (- 4 * p.x);

//    der [2].x = LAMBDA_EVAL (-1 + 4 * p.x);
//    der [2].y = LAMBDA_EVAL (0);

//    der [3].x = LAMBDA_EVAL (4. * p.y);
//    der [3].y = LAMBDA_EVAL (4. * p.x);

//    der [4].x = LAMBDA_EVAL (0);
//    der [4].y = LAMBDA_EVAL (-1 + 4 * p.y);

//    der [5].x = LAMBDA_EVAL (- 4 * p.y);
//    der [5].y = LAMBDA_EVAL (4 * ((1. - p.x - p.y) - p.y));

    return;
}
