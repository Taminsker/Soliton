#include "functions.h"

FunAnalytic      fun_analytic;
FunAnalytic_1    fun_analytic_1;
FunAnalytic_2    fun_analytic_2;

FunNeumann       fun_neumann;
FunDirichlet     fun_dirichlet;
FunSecondMember  fun_secondmember;
FunJumpZeta      fun_jumpzeta;

real_t           coeffA = 0.0;
real_t           coeffT = 0.0;
real_t           coeffLambda = 0.0;

real_t
FunAnalytic::ToReal (Point *p, real_t t, Cell *)
{
    real_t omega = 2. * M_PI / coeffT;
    real_t c = coeffLambda / coeffT;
    real_t k = omega / c;

//    return coeffA * std::sin(omega * t) * std::sin(k * p->x);

    return coeffA * std::cos (omega * t - k * p->x);
}

real_t
FunAnalytic_1::ToReal (Point *p, real_t t, Cell *)
{
    real_t omega = 2. * M_PI / coeffT;
    real_t c = coeffLambda / coeffT;
    real_t k = omega / c;

//    return coeffA * std::sin(omega * t) * std::sin(k * p->x);

    return - coeffA * omega * std::sin (omega * t - k * p->x);
}

real_t
FunAnalytic_2::ToReal (Point *p, real_t t, Cell *cell)
{
    real_t omega = 2. * M_PI / coeffT;
//    real_t c = coeffLambda / coeffT;
//    real_t k = omega / c;

//    return coeffA * std::sin(omega * t) * std::sin(k * p->x);

    return omega * omega * fun_analytic.ToReal (p, t, cell);
}

Point
FunNeumann::To3DPoint (Point *, real_t, Cell *)
{
    return Point();
}

real_t
FunDirichlet::ToReal (Point *p, real_t t, Cell * c)
{
    return fun_analytic.ToReal (p, t, c);
}

real_t
FunSecondMember::ToReal (Point * p, real_t t, Cell *cel)
{
//    return 2. * (p->x - p->x * p->x + p->y - p->y * p->y);
//        return 0.;
    real_t omega = 2. * M_PI / coeffT;
    real_t c = coeffLambda / coeffT;
    real_t k = omega / c;
    return k*k * fun_analytic.ToReal (p, t, cel);

    //  return 0.;
    //  return 1.;
    //    return 0;
    //    return 2. * (p.x - p.x * p.x + p.y - p.y * p.y);
    //  return -2 * (1. + 1.);
    //    return 2. * M_PI * M_PI * fun_analytic (p, t);
}

real_t
FunJumpZeta::ToReal (Point *, real_t, Cell *)
{
    return 0.0;
}
