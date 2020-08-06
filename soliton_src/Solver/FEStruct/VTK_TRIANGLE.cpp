#include "festruct.h"

template <>
void LagrangeFE<VTK_TRIANGLE>::Build ()
{
    npts = 3;
    coor = {{0, 0}, {1, 0}, {0, 1}};
    volRef = 0.5;

    fun.resize (npts);
    der.resize (npts);

    fun [0]   = LAMBDA_EVAL (1. - p.x - p.y);
    fun [1]   = LAMBDA_EVAL (p.x);
    fun [2]   = LAMBDA_EVAL (p.y);

    der [0].x = LAMBDA_EVAL (-1.);
    der [0].y = LAMBDA_EVAL (-1.);

    der [1].x = LAMBDA_EVAL (1.);
    der [1].y = LAMBDA_EVAL (0.);

    der [2].x = LAMBDA_EVAL (0.);
    der [2].y = LAMBDA_EVAL (1.);


    return;
}


template <>
std::matrix<double> LagrangeFE <VTK_TRIANGLE>::Jac (Cell* cell)
{
    auto P = cell->GetPoints ();

    if (P->size () != 3)
        return Math::Zeros (3);

    std::matrix<double> mat = Math::Ones (3);

    mat[0][0] = P->at(1)->x - P->at(0)->x;
    mat[1][0] = P->at(1)->y - P->at(0)->y;
    mat[2][0] = P->at(1)->z - P->at(0)->z;

    mat[0][1] = P->at(2)->x - P->at(0)->x;
    mat[1][1] = P->at(2)->y - P->at(0)->y;
    mat[2][1] = P->at(2)->z - P->at(0)->z;

    return mat;
}

template <>
Point LagrangeFE <VTK_TRIANGLE>::B (Cell* cell)
{
    if (cell->GetPoints ()->size () != 1)
        return Point();

    return *cell->GetPoints ()->at (0);
}
