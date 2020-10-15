#include "functions.h"

double fun_analytic (Point p, double t)
{
    VOID_USE(t);
    return (p.x - p.x * p.x) * (p.y - p.y * p.y);
//        return (p.x * p.x - p.y * p.y);
//        return p.x + p.y;

}

double fun_neumann (Point p, double t)
{
    VOID_USE(p);
    VOID_USE(t);

    return (p.x - p.x * p.x) * (1. - 2. * p.y);
    //    return (p.x * p.x + p.y * p.y);
}

double fun_dirichlet (Point p, double t)
{
    VOID_USE(p);
    VOID_USE(t);

    return fun_analytic (p, t);
}

double fun_secondmember (Point p, double t)
{
    VOID_USE(t);
    VOID_USE(p);

    return 2. * (p.x - p.x * p.x + p.y - p.y * p.y);
//        return 0.;
//    return 0.;
}
