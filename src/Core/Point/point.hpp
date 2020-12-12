#ifndef SRC_CORE_POINT_POINT_HPP
#define SRC_CORE_POINT_POINT_HPP

#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "../../Algorithms/Math/math.hpp"
#include "../../solitonheader.hpp"

#define NONE_ID_SELECTED 0

class Cell;

class Point
{
public:
    real_t x = 0x0;
    real_t y = 0x0;
    real_t z = 0x0;

    Point ();
    Point (real_t a, real_t b = 0., real_t c = 0.);
    Point (const Point & p);
    ~Point ();

    SOLITON_INLINE
    Point &
    operator= (std::initializer_list<real_t> ilist)
    {
        auto iter = ilist.begin ();
        x         = (ilist.size () >= 1 ? *iter : 0.);
        iter++;
        y = (ilist.size () >= 2 ? *iter : 0.);
        iter++;
        z = (ilist.size () >= 3 ? *iter : 0.);

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator= (std::initializer_list<int> ilist)
    {
        auto iter = ilist.begin ();
        x         = real_t (ilist.size () >= 1 ? *iter : 0.);
        iter++;
        y = real_t (ilist.size () >= 2 ? *iter : 0.);
        iter++;
        z = real_t (ilist.size () >= 3 ? *iter : 0.);

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator= (const Point & p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
        //  m_cells = p.m_cells;
        //  m_listNeighbours = p.m_listNeighbours;

        //  if (p.m_globalIndex != NONE_ID_SELECTED)
        //    m_globalIndex = p.m_globalIndex;

        return *this;
    }

    SOLITON_INLINE
    Point *
    SetGlobalIndex (ul_t index)
    {
        m_globalIndex = index;
        return this;
    }

    SOLITON_INLINE
    ul_t
    GetGlobalIndex () const
    {
        return m_globalIndex;
    }

    SOLITON_INLINE
    void
    ClearListNeighbours ()
    {
        m_listNeighbours.clear ();
        return;
    }

    void AddPointNeighbour (Point * p);

    SOLITON_INLINE
    std::vector<Point *>
    GetListNeighbours ()
    {
        return m_listNeighbours;
    }

    void SetListNeighbours (std::vector<Point *> & list);
    void RemoveThisNeighbourPoint (Point * p);
    void LinkToCell (Cell * c);
    void UnlinkToCell (Cell * c);

    SOLITON_INLINE
    std::vector<Cell *>
    GetLinkedCell () const
    {
        return m_cells;
    }

    Point * DetachFromAll ();

    SOLITON_INLINE
    Point &
    operator+= (Point && p)
    {
        x += p.x;
        y += p.y;
        z += p.z;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator-= (Point && p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator*= (Point && p)
    {
        x *= p.x;
        y *= p.y;
        z *= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator+= (Point & p)
    {
        x += p.x;
        y += p.y;
        z += p.z;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator-= (Point & p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator*= (Point & p)
    {
        x *= p.x;
        y *= p.y;
        z *= p.z;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator+= (real_t value)
    {
        x += value;
        y += value;
        z += value;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator-= (real_t value)
    {
        x -= value;
        y -= value;
        z -= value;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator*= (real_t value)
    {
        x *= value;
        y *= value;
        z *= value;

        return *this;
    }

    SOLITON_INLINE
    Point &
    operator/= (real_t value)
    {
        if (std::abs (value) < EPSILON)
            return *this;

        x /= value;
        y /= value;
        z /= value;

        return *this;
    }

    SOLITON_INLINE
    std::vector<real_t>
    data ()
    {
        return {x, y, z};
    }

    SOLITON_INLINE
    real_t
    EuclidianNorm ()
    {
        return std::sqrt (x * x + y * y + z * z);
    }

    SOLITON_INLINE
    Point &
    Normalize ()
    {
        real_t norm = this->EuclidianNorm ();
        if (std::abs (norm) < EPSILON)
            return *this;

        this->x /= norm;
        this->y /= norm;
        this->z /= norm;

        return *this;
    }

    friend std::ostream & operator<< (std::ostream & out, const Point & p);
    friend bool           CompareIdx (const Point & a, const Point & b);

    friend real_t    EuclidianDist (const Point & a, const Point & b);
    friend Point     CrossProduct (const Point & a, const Point & b);
    friend Matrix3x3 OuterProduct (const Point & a, const Point & b);
    friend Point     operator* (real_t value, const Point & p);
    friend Point     operator* (const Point & p, real_t value);
    friend Point     operator* (const Point & a, const Point & b);
    friend Point     operator* (const Matrix3x3 & mat, const Point & p);

    friend Point operator+ (const Point & a, const Point & b);
    friend Point operator- (const Point & a, const Point & b);

    friend real_t operator| (const Point & a, const Point & b);

    friend Point operator+ (const Point & p, real_t value);
    friend Point operator+ (real_t value, const Point & p);

    friend Point operator- (const Point & p, real_t value);
    friend Point operator- (real_t value, const Point & p);

    friend Point operator/ (const Point & p, real_t value);
    friend Point operator/ (real_t value, const Point & p);

    friend bool operator< (const Point & a, const Point & b);
    friend bool operator> (const Point & a, const Point & b);

    friend bool operator<= (const Point & a, const Point & b);
    friend bool operator>= (const Point & a, const Point & b);

    friend bool operator== (const Point & a, const Point & b);
    friend bool operator!= (const Point & a, const Point & b);

protected:
    std::vector<Cell *>  m_cells;
    std::vector<Point *> m_listNeighbours;
    ul_t                 m_globalIndex;
};

std::ostream & operator<< (std::ostream & out, const Point & p);
std::ostream & operator<< (std::ostream & out, const std::vector<Point *> vec);

bool CompareIdx (const Point & a, const Point & b);

Point operator* (real_t value, const Point & p);
Point operator* (const Point & p, real_t value);
Point operator* (const Point & a, const Point & b);
Point operator* (const Matrix3x3 & mat, const Point & p);

Point  operator+ (const Point & a, const Point & b);
Point  operator- (const Point & a, const Point & b);
real_t operator| (const Point & a, const Point & b);
Point  operator+ (const Point & p, real_t value);
Point  operator+ (real_t value, const Point & p);
Point  operator- (const Point & p, real_t value);
Point  operator- (real_t value, const Point & p);
Point  operator/ (const Point & p, real_t value);
Point  operator/ (real_t value, const Point & p);
bool   operator< (const Point & a, const Point & b);
bool   operator> (const Point & a, const Point & b);
bool   operator<= (const Point & a, const Point & b);
bool   operator>= (const Point & a, const Point & b);
bool   operator== (const Point & a, const Point & b);
bool   operator!= (const Point & a, const Point & b);

real_t    EuclidianDist (const Point & a, const Point & b);
Point     CrossProduct (const Point & a, const Point & b);
Matrix3x3 OuterProduct (const Point & a, const Point & b);

void InitPointVector (std::vector<Point *> * vec, ul_t size, Point && data = Point (0, 0, 0));
#endif /* SRC_CORE_POINT_POINT_HPP */
