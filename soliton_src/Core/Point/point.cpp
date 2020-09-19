#include "point.h"
#include "../Cell/cell.h"

Point::Point () :
    x (0.), y (0.), z (0.),
    m_cells ({}),
    m_listNeighbours ({}),
    m_globalIndex (NONE_ID_SELECTED) // Il n'a pas encore d'indice global
{}

Point::Point (const Point &p) :
    x (p.x), y (p.y), z (p.z),
    m_cells (p.m_cells),
    m_listNeighbours (p.m_listNeighbours),
    m_globalIndex (p.m_globalIndex)
{}

Point::Point (double a, double b, double c) :
    x (a), y (b), z (c),
    m_cells ({}),
    m_listNeighbours ({}),
    m_globalIndex (NONE_ID_SELECTED) // Il n'a pas encore d'indice global
{}

Point& Point::operator= (const Point& p)
{
    x = p.x;
    y = p.y;
    z = p.z;
    m_cells  = p.m_cells;
    m_listNeighbours  = p.m_listNeighbours;

    if (p.m_globalIndex != NONE_ID_SELECTED)
        m_globalIndex = p.m_globalIndex;

    return *this;
}

Point& Point::operator= (std::initializer_list<double> ilist)
{
    auto iter = ilist.begin ();
    x = (ilist.size () >= 1 ? *iter : 0.);
    iter++;
    y = (ilist.size () >= 2 ? *iter : 0.);
    iter++;
    z = (ilist.size () >= 3 ? *iter : 0.);

    return *this;
}

Point& Point::operator= (std::initializer_list<int> ilist)
{
    auto iter = ilist.begin ();
    x = double (ilist.size () >= 1 ? *iter : 0.);
    iter++;
    y = double (ilist.size () >= 2 ? *iter : 0.);
    iter++;
    z = double (ilist.size () >= 3 ? *iter : 0.);

    return *this;
}

Point::~Point ()
{}

Point* Point::SetGlobalIndex (int index)
{
    m_globalIndex = index;
    return this;
}

int Point::GetGlobalIndex () const
{
    return m_globalIndex;
}

void Point::AddPointNeighbour (Point* p)
{
    for (auto up : m_listNeighbours)
    {
        if (*up == *p)
            return;
    }

    m_listNeighbours.push_back (p);

    return;
}

void Point::RemoveThisNeighbourPoint (Point* p)
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

void Point::ClearListNeighbours ()
{
    m_listNeighbours.clear ();
    return;
}


std::vector<Point*> Point::GetListNeighbours ()
{
    return m_listNeighbours;
}

void Point::SetListNeighbours (std::vector <Point *>& list)
{
    for (Point * p : list)
        if (p == nullptr)
            return;

    m_listNeighbours = list;

    return;
}

void Point::LinkToCell (Cell *c)
{
    for (Cell* tc : m_cells)
    {
        if (tc == c)
            return;
    }

    m_cells.push_back (c);

    return;
}

void Point::UnlinkToCell (Cell *c)
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

std::vector <Cell*> Point::GetLinkedCell () const
{
    return m_cells;
}

Point* Point::DetachFromAll ()
{
    for (Cell* c : m_cells)
        c->RemovePoint (this);

    for (Point* p : m_listNeighbours)
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

Point& Point::operator+= (Point&& p)
{
    x += p.x;
    y += p.y;
    z += p.z;

    return *this;
}

Point& Point::operator-= (Point&& p)
{
    x -= p.x;
    y -= p.y;
    z -= p.z;

    return *this;
}

Point& Point::operator*= (Point&& p)
{
    x *= p.x;
    y *= p.y;
    z *= p.z;

    return *this;
}

Point& Point::operator+= (Point& p)
{
    x += p.x;
    y += p.y;
    z += p.z;

    return *this;
}

Point& Point::operator-= (Point& p)
{
    x -= p.x;
    y -= p.y;
    z -= p.z;

    return *this;
}

Point& Point::operator*= (Point& p)
{
    x *= p.x;
    y *= p.y;
    z *= p.z;

    return *this;
}

Point& Point::operator+= (double value)
{
    x += value;
    y += value;
    z += value;

    return *this;
}

Point& Point::operator-= (double value)
{
    x -= value;
    y -= value;
    z -= value;

    return *this;
}

Point& Point::operator*= (double value)
{
    x *= value;
    y *= value;
    z *= value;

    return *this;
}

Point& Point::operator/= (double value)
{
    if (std::abs (value) < 1e-20)
        return *this;

    x /= value;
    y /= value;
    z /= value;

    return *this;
}

std::vector<double> Point::data ()
{
    return {x, y, z};
}

Point& Point::Normalize ()
{
    double norm = this->EuclidianNorm ();
    if (std::abs (norm) < 1e-15)
        return *this;

    this->x /= norm;
    this->y /= norm;
    this->z /= norm;

    return *this;
}

double Point::EuclidianNorm ()
{
    return std::sqrt (x * x + y * y + z * z);
}


std::ostream & operator<< (std::ostream &out, const Point &p)
{
    out << SPACE << p.x << " "
        << SPACE << p.y << " "
        << SPACE << p.z << " ";
    return out;
}

std::ostream & operator<< (std::ostream &out, const std::vector<Point*> vec)
{
    for (std::size_t i = 0; i < vec.size (); ++i)
        out << *(vec.at (i)) << std::endl;
    return out;
}

Point operator* (double value, const Point& p)
{
    return {value * p.x, value * p.y, value * p.z};
}

Point operator* (const Point& p, double value)
{
    return value * p;
}

Point operator* (const Point& a, const Point& b)
{
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}

Point operator* (const Matrix3x3& A, const Point& p)
{
    Point out;
    out.x = A.coeffRef (0, 0) * p.x + A.coeffRef (0, 1) * p.y + A.coeffRef (0, 2) * p.z;
    out.y = A.coeffRef (1, 0) * p.x + A.coeffRef (1, 1) * p.y + A.coeffRef (1, 2) * p.z;
    out.z = A.coeffRef (2, 0) * p.x + A.coeffRef (2, 1) * p.y + A.coeffRef (2, 2) * p.z;

    return out;
}
Point operator+ (const Point& a, const Point& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Point operator- (const Point& a, const Point& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

double operator| (const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point operator+ (double value, const Point& p)
{
    return {value + p.x, value + p.y, value + p.z};
}

Point operator+ (const Point& p, double value)
{
    return value + p;
}

Point operator- (double value, const Point& p)
{
    return {value - p.x, value - p.y, value - p.z};
}

Point operator- (const Point& p, double value)
{
    return (-1 * value) + p;
}

Point operator/ (double value, const Point& p)
{
    return {value / p.x, value / p.y, value / p.z};
}

Point operator/ (const Point& p, double value)
{
    return (1. / value) * p;
}

bool operator< (const Point& a, const Point& b)
{
    return ((a.x < b.x) && (a.y < b.y) && (a.z < b.z));
}

bool operator> (const Point& a, const Point& b)
{
    return ((a.x > b.x) && (a.y > b.y) && (a.z > b.z));
}

bool operator<= (const Point& a, const Point& b)
{
    return ((a.x <= b.x) && (a.y <= b.y) && (a.z <= b.z));
}

bool operator>= (const Point& a, const Point& b)
{
    return ((a.x >= b.x) && (a.y >= b.y) && (a.z >= b.z));
}

bool operator== (const Point& a, const Point& b)
{
    bool bx = (fabs(a.x - b.x) < EPSILON);
    bool by = (fabs(a.y - b.y) < EPSILON);
    bool bz = (fabs(a.z - b.z) < EPSILON);

    return (bx && by && bz);
}

bool operator!= (const Point& a, const Point& b)
{
    bool bx = (fabs(a.x - b.x) > EPSILON);
    bool by = (fabs(a.y - b.y) > EPSILON);
    bool bz = (fabs(a.z - b.z) > EPSILON);

    return (bx && by && bz);
}


double EuclidianDist (const Point& a, const Point& b)
{
    double x = a.x - b.x;
    double y = a.y - b.y;
    double z = a.z - b.z;

    return std::sqrt (x * x + y * y + z * z);
}

Point CrossProduct (const Point& a, const Point& b)
{
    Point out;
    out.x = a.y * b.z - a.z * b.y;
    out.y = a.z * b.x - a.x * b.z;
    out.z = a.x * b.y - a.y * b.x;
    return out;
}

Matrix3x3 OuterProduct (const Point& a, const Point& b)
{
    Matrix3x3 out;
    out.setZero ();

    out.row (0) << a.x * b.x, a.x * b.y, a.x * b.z;
    out.row (1) << a.y * b.x, a.y * b.y, a.y * b.z;
    out.row (2) << a.z * b.x, a.z * b.y, a.z * b.z;

    return out;
}

void InitPointVector (std::vector<Point*>* vec, std::size_t size, Point&& data)
{
    vec->resize (size);
    for (std::size_t i = 0; i < size; ++i)
    {
        vec->at (i) = new Point();
        *vec->at (i) = data;
    }

    return;
}



























