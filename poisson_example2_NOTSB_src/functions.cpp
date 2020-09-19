#include "functions.h"

double u (Point p, double t)
{
    (void)t;
    return std::exp (p.x + p.y) * std::sin (p.x) * std::cos (p.y);
}

double g (Point p, double t)
{
    (void)p;
    (void)t;

    return u (p, t);
}

double h (Point p, double t)
{
    (void)p;
    (void)t;

    return u (p, t);

}

double f (Point p, double t)
{
    (void)t;
    (void)p;

    return -2. * std::exp (p.x + p.y) * std::cos (p.x + p.y);
}
