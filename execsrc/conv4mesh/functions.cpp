#include "functions.h"

FunAnalytic      fun_analytic;
FunNeumann       fun_neumann;
FunDirichlet     fun_dirichlet;
FunSecondMember  fun_secondmember;

real_t
FunAnalytic::ToReal (Point *p, real_t, Cell *)
{
    return (p->x - p->x * p->x) * (p->y - p->y * p->y);
    //  return 1. * p.x * p.x + 1. * p.y * p.y + 1. * p.x * p.y + 1. * p.x + 1. * p.y + 1.;
    //    return p.x + p.y;
    //  return -0.25 * p.x * p.x - 0.25 * p.y * p.y;
    //  return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
    //    return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
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
FunSecondMember::ToReal (Point *p, real_t, Cell *)
{
    return 2. * (p->x - p->x * p->x + p->y - p->y * p->y);
    //    return 0.;
    //  return 0.;
    //  return 1.;
    //    return 0;
//    return 2. * (p.x - p.x * p.x + p.y - p.y * p.y);
    //  return -2 * (1. + 1.);
    //    return 2. * M_PI * M_PI * fun_analytic (p, t);
}
