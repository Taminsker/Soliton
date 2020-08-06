#include "festruct.h"

template <>
void LagrangeFE<VTK_VERTEX>::Build ()
{
    npts = 1;
    coor = {{0}};
    volRef = 0;

    fun.resize (npts);
    der.resize (npts);

    fun [0] = LAMBDA_EVAL (0.);

    der [0] = Functor_L2P();

    return;
}

template <>
std::matrix<double> LagrangeFE<VTK_VERTEX>::Jac (Cell* cell)
{
    if (cell->GetPoints ()->size () != 1)
        return Math::Zeros (3);
    return Math::Ones (3);
}

template <>
Point LagrangeFE<VTK_VERTEX>::B (Cell* cell)
{
    if (cell->GetPoints ()->size () != 1)
        return Point();

    return *cell->GetPoints ()->at (0);
}
