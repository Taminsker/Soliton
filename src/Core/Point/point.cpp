#include "point.hpp"

#include "../../Algorithms/Math/math.hpp"
#include "../Cell/cell.hpp"

Point::Point () : x (0.),
                  y (0.),
                  z (0.),
                  m_cells ({}),
                  m_listNeighbours ({}),
                  m_globalIndex (NONE_ID_SELECTED)  // Il n'a pas encore d'indice global
{
}

Point::Point (const Point & p) : x (p.x),
                                 y (p.y),
                                 z (p.z),
                                 m_cells (p.m_cells),
                                 m_listNeighbours (p.m_listNeighbours),
                                 m_globalIndex (p.m_globalIndex)
{
}

Point::Point (real_t a, real_t b, real_t c) : x (a),
                                              y (b),
                                              z (c),
                                              m_cells ({}),
                                              m_listNeighbours ({}),
                                              m_globalIndex (NONE_ID_SELECTED)  // Il n'a pas encore d'indice global
{
}

Point::~Point ()
{
}

void
Point::AddPointNeighbour (Point * p)
{
    for (auto up : m_listNeighbours)
    {
        if (*up == *p)
            return;
    }

    m_listNeighbours.push_back (p);

    return;
}

void
Point::RemoveThisNeighbourPoint (Point * p)
{
    auto it = m_listNeighbours.begin ();

    while (it != m_listNeighbours.end ())
    {
        if (p == *it)
            it = m_listNeighbours.erase (it);
        else
            it++;
    }

    return;
}

void
Point::SetListNeighbours (std::vector<Point *> & list)
{
    for (Point * p : list)
        if (p == nullptr)
            return;

    m_listNeighbours = list;

    return;
}

void
Point::LinkToCell (Cell * c)
{
    for (Cell * tc : m_cells)
    {
        if (tc == c)
            return;
    }

    m_cells.push_back (c);

    return;
}

void
Point::UnlinkToCell (Cell * c)
{
    auto it = m_cells.begin ();
    while (it != m_cells.end ())
    {
        if (*it == c)
            m_cells.erase (it);
        else
            it++;
    }

    return;
}

Point *
Point::DetachFromAll ()
{
    for (Cell * c : m_cells)
        c->RemovePoint (this);

    for (Point * p : m_listNeighbours)
        p->RemoveThisNeighbourPoint (this);
    if (m_listNeighbours.size () == 2)
    {
        m_listNeighbours.at (0)->AddPointNeighbour (m_listNeighbours.at (1));
        m_listNeighbours.at (1)->AddPointNeighbour (m_listNeighbours.at (0));
    }

    m_listNeighbours.clear ();
    m_cells.clear ();

    return this;
}

std::ostream &
operator<< (std::ostream & out, const Point & p)
{
    out << SPACE << p.x << " "
        << SPACE << p.y << " "
        << SPACE << p.z << " ";
    return out;
}

std::ostream &
operator<< (std::ostream & out, const std::vector<Point *> vec)
{
    for (ul_t i = 0; i < vec.size (); ++i)
        out << *(vec.at (i)) << std::endl;
    return out;
}

bool
CompareIdx (const Point & a, const Point & b)
{
    return a.GetGlobalIndex () > b.GetGlobalIndex ();
}

Point
operator* (real_t value, const Point & p)
{
    return {value * p.x, value * p.y, value * p.z};
}

Point
operator* (const Point & p, real_t value)
{
    return value * p;
}

Point
operator* (const Point & a, const Point & b)
{
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}

Point
operator* (const Matrix3x3 & A, const Point & p)
{
    Point out;
    out.x = A.coeffRef (0, 0) * p.x + A.coeffRef (0, 1) * p.y + A.coeffRef (0, 2) * p.z;
    out.y = A.coeffRef (1, 0) * p.x + A.coeffRef (1, 1) * p.y + A.coeffRef (1, 2) * p.z;
    out.z = A.coeffRef (2, 0) * p.x + A.coeffRef (2, 1) * p.y + A.coeffRef (2, 2) * p.z;

    return out;
}
Point
operator+ (const Point & a, const Point & b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Point
operator- (const Point & a, const Point & b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

real_t
operator| (const Point & a, const Point & b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point
operator+ (real_t value, const Point & p)
{
    return {value + p.x, value + p.y, value + p.z};
}

Point
operator+ (const Point & p, real_t value)
{
    return value + p;
}

Point
operator- (real_t value, const Point & p)
{
    return {value - p.x, value - p.y, value - p.z};
}

Point
operator- (const Point & p, real_t value)
{
    return (-1 * value) + p;
}

Point
operator/ (real_t value, const Point & p)
{
    return {value / p.x, value / p.y, value / p.z};
}

Point
operator/ (const Point & p, real_t value)
{
    return (1. / value) * p;
}

bool
operator< (const Point & a, const Point & b)
{
    return ((a.x < b.x) && (a.y < b.y) && (a.z < b.z));
}

bool
operator> (const Point & a, const Point & b)
{
    return ((a.x > b.x) && (a.y > b.y) && (a.z > b.z));
}

bool
operator<= (const Point & a, const Point & b)
{
    return ((a.x <= b.x) && (a.y <= b.y) && (a.z <= b.z));
}

bool
operator>= (const Point & a, const Point & b)
{
    return ((a.x >= b.x) && (a.y >= b.y) && (a.z >= b.z));
}

bool
operator== (const Point & a, const Point & b)
{
    bool bx = (fabs (a.x - b.x) < EPSILON);
    bool by = (fabs (a.y - b.y) < EPSILON);
    bool bz = (fabs (a.z - b.z) < EPSILON);

    return (bx && by && bz);
}

bool
operator!= (const Point & a, const Point & b)
{
    bool bx = (fabs (a.x - b.x) > EPSILON);
    bool by = (fabs (a.y - b.y) > EPSILON);
    bool bz = (fabs (a.z - b.z) > EPSILON);

    return (bx && by && bz);
}

real_t
EuclidianDist (const Point & a, const Point & b)
{
    real_t x = a.x - b.x;
    real_t y = a.y - b.y;
    real_t z = a.z - b.z;

    return std::sqrt (x * x + y * y + z * z);
}

Point
CrossProduct (const Point & a, const Point & b)
{
    Point out;
    out.x = a.y * b.z - a.z * b.y;
    out.y = a.z * b.x - a.x * b.z;
    out.z = a.x * b.y - a.y * b.x;
    return out;
}

Matrix3x3
OuterProduct (const Point & a, const Point & b)
{
    Matrix3x3 out;
    out.setZero ();

    out.row (0) << a.x * b.x, a.x * b.y, a.x * b.z;
    out.row (1) << a.y * b.x, a.y * b.y, a.y * b.z;
    out.row (2) << a.z * b.x, a.z * b.y, a.z * b.z;

    return out;
}

void
InitPointVector (std::vector<Point *> * vec, ul_t size, Point && data)
{
    vec->resize (size);
    for (ul_t i = 0; i < size; ++i)
    {
        vec->at (i)  = new Point ();
        *vec->at (i) = data;
    }

    return;
}
