
#include "festruct.h"

template <>
void LagrangeFE<VTK_LINE>::Build ()
{
    npts = 2;
    coor = {{-1.}, {1.}};
    volRef = 2.;

    fun.resize (npts);
    der.resize (npts);

    fun [0]   = LAMBDA_EVAL (2. * (1. - p.x));
    fun [1]   = LAMBDA_EVAL (2. * (1. + p.x));

    der [0].x = LAMBDA_EVAL (-2.);
    der [1].x = LAMBDA_EVAL (2.);

    return;
}

template <>
std::matrix<double> LagrangeFE<VTK_LINE>::Jac (Cell* cell)
{
    auto P = cell->GetPoints ();

    if (P->size () != 2)
        return Math::Zeros (3);

    std::matrix<double> mat = Math::Ones (3);

    mat[0][0] = 0.5 * (P->at(1)->x - P->at(0)->x);
    mat[1][0] = 0.5 * (P->at(1)->y - P->at(0)->y);
    mat[2][0] = 0.5 * (P->at(1)->z - P->at(0)->z);

    return mat;
}

template <>
Point LagrangeFE<VTK_LINE>::B (Cell* cell)
{
    auto P = cell->GetPoints ();

    if (P->size () != 1)
        return Point();

    return 0.5 * (*P->at(0) + *P->at(1));
}
