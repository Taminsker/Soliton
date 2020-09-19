#include "functions.h"

double fun_analytic (Point p, double t)
{
    (void)t;
    return (p.x - p.x * p.x) * (p.y - p.y * p.y);
//    return p.x + p.y;
}

double fun_neumann (Point p, double t)
{
    (void)p;
    (void)t;

    return 0.;
//        return (p.x * p.x + p.y * p.y);
}

double fun_dirichlet (Point p, double t)
{
    (void)p;
    (void)t;

    return fun_analytic (p, t);
}

double fun_secondmember (Point p, double t)
{
    (void)t;
    (void)p;

//    return 0;
    return 2. * (p.x - p.x * p.x + p.y - p.y * p.y);
}
