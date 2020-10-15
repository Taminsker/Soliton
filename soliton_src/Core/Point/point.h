#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <initializer_list>
#include <utility>

#include <Algorithms/Math/math.h>
#include <ProgDef/proddef.h>

#define NONE_ID_SELECTED -2

class Cell;

class Point
{
public:
    double x = 0x0;
    double y = 0x0;
    double z = 0x0;

    Point ();
    Point (double a, double b = 0., double c = 0.);
    Point (const Point &p);
    ~Point ();

    SOLITON_INLINE
    Point& operator= (std::initializer_list<double> ilist)
    {
        auto iter = ilist.begin ();
        x = (ilist.size () >= 1 ? *iter : 0.);
        iter++;
        y = (ilist.size () >= 2 ? *iter : 0.);
        iter++;
        z = (ilist.size () >= 3 ? *iter : 0.);

        return *this;
    }

    SOLITON_INLINE
    Point& operator= (std::initializer_list<int> ilist)
    {
        auto iter = ilist.begin ();
        x = double (ilist.size () >= 1 ? *iter : 0.);
        iter++;
        y = double (ilist.size () >= 2 ? *iter : 0.);
        iter++;
        z = double (ilist.size () >= 3 ? *iter : 0.);

        return *this;
    }

    SOLITON_INLINE
    Point& operator= (const Point& p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
        //    m_cells  = p.m_cells;
        //    m_listNeighbours  = p.m_listNeighbours;

        //    if (p.m_globalIndex != NONE_ID_SELECTED)
        //        m_globalIndex = p.m_globalIndex;

        return *this;
    }

    SOLITON_INLINE
    Point* SetGlobalIndex (int index)
    {
        m_globalIndex = index;
        return this;
    }

    SOLITON_INLINE
    int GetGlobalIndex () const
    {
        return m_globalIndex;
    }

    SOLITON_INLINE
    void ClearListNeighbours ()
    {
        m_listNeighbours.clear ();
        return;
    }

    void AddPointNeighbour (Point* p);

    SOLITON_INLINE
    std::vector <Point *> GetListNeighbours ()
    {
        return m_listNeighbours;
    }

    void SetListNeighbours (std::vector <Point *>& list);
    void RemoveThisNeighbourPoint (Point* p);
    void LinkToCell (Cell* c);
    void UnlinkToCell (Cell* c);

    SOLITON_INLINE
    std::vector <Cell*> GetLinkedCell () const
    {
        return m_cells;
    }

    Point*DetachFromAll ();

    SOLITON_INLINE
    Point& operator+= (Point&& p)
    {
        x += p.x;
        y += p.y;
        z += p.z;

        return *this;
    }

    SOLITON_INLINE
    Point& operator-= (Point&& p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point& operator*= (Point&& p)
    {
        x *= p.x;
        y *= p.y;
        z *= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point& operator+= (Point& p)
    {
        x += p.x;
        y += p.y;
        z += p.z;

        return *this;
    }

    SOLITON_INLINE
    Point& operator-= (Point& p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point& operator*= (Point& p)
    {
        x *= p.x;
        y *= p.y;
        z *= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point& operator+= (double value)
    {
        x += value;
        y += value;
        z += value;

        return *this;
    }

    SOLITON_INLINE
    Point& operator-= (double value)
    {
        x -= value;
        y -= value;
        z -= value;

        return *this;
    }

    SOLITON_INLINE
    Point& operator*= (double value)
    {
        x *= value;
        y *= value;
        z *= value;

        return *this;
    }

    SOLITON_INLINE
    Point& operator/= (double value)
    {
        if (std::abs (value) < EPSILON)
            return *this;

        x /= value;
        y /= value;
        z /= value;

        return *this;
    }

    SOLITON_INLINE
    std::vector<double> data ()
    {
        return {x, y, z};
    }

    SOLITON_INLINE
    double EuclidianNorm ()
    {
        return std::sqrt (x * x + y * y + z * z);
    }

    SOLITON_INLINE
    Point& Normalize ()
    {
        double norm = this->EuclidianNorm ();
        if (std::abs (norm) < EPSILON)
            return *this;

        this->x /= norm;
        this->y /= norm;
        this->z /= norm;

        return *this;
    }

    friend std::ostream & operator<< (std::ostream &out, const Point &p);
    friend bool CompareIdx (const Point& a, const Point& b);

    friend double EuclidianDist (const Point& a, const Point& b);
    friend Point CrossProduct (const Point& a, const Point& b);
    friend Matrix3x3 OuterProduct (const Point& a, const Point& b);
    friend Point operator* (double value, const Point& p);
    friend Point operator* (const Point& p, double value);
    friend Point operator* (const Point& a, const Point& b);
    friend Point operator* (const Matrix3x3& mat, const Point& p);

    friend Point operator+ (const Point& a, const Point& b);
    friend Point operator- (const Point& a, const Point& b);

    friend double operator| (const Point& a, const Point& b);

    friend Point operator+ (const Point& p, double value);
    friend Point operator+ (double value, const Point& p);

    friend Point operator- (const Point& p, double value);
    friend Point operator- (double value, const Point& p);

    friend Point operator/ (const Point& p, double value);
    friend Point operator/ (double value, const Point& p);

    friend bool operator< (const Point& a, const Point& b);
    friend bool operator> (const Point& a, const Point& b);

    friend bool operator<= (const Point& a, const Point& b);
    friend bool operator>= (const Point& a, const Point& b);

    friend bool operator== (const Point& a, const Point& b);
    friend bool operator!= (const Point& a, const Point& b);

protected:
    std::vector<Cell *> m_cells;
    std::vector<Point *> m_listNeighbours;
    int m_globalIndex;
};

std::ostream& operator<< (std::ostream &out, const Point &p);
std::ostream & operator<< (std::ostream &out, const std::vector<Point*> vec);

bool CompareIdx (const Point& a, const Point& b);

Point operator* (double value, const Point& p);
Point operator* (const Point& p, double value);
Point operator* (const Point& a, const Point& b);
Point operator* (const Matrix3x3& mat, const Point& p);

Point operator+ (const Point& a, const Point& b);
Point operator- (const Point& a, const Point& b);
double operator| (const Point& a, const Point& b);
Point operator+ (const Point& p, double value);
Point operator+ (double value, const Point& p);
Point operator- (const Point& p, double value);
Point operator- (double value, const Point& p);
Point operator/ (const Point& p, double value);
Point operator/ (double value, const Point& p);
bool operator< (const Point& a, const Point& b);
bool operator> (const Point& a, const Point& b);
bool operator<= (const Point& a, const Point& b);
bool operator>= (const Point& a, const Point& b);
bool operator== (const Point& a, const Point& b);
bool operator!= (const Point& a, const Point& b);

double EuclidianDist (const Point& a, const Point& b);
Point CrossProduct (const Point& a, const Point& b);
Matrix3x3 OuterProduct (const Point& a, const Point& b);

void InitPointVector (std::vector<Point*>* vec, std::size_t size, Point&& data = Point(0, 0, 0));
#endif // POINT_H
