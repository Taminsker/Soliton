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
    double x = 0.;
    double y = 0.;
    double z = 0.;

    Point ();
    Point (const Point &p);
    Point (double a, double b = 0., double c = 0.);
    Point& operator= (const Point& p);
    Point& operator= (std::initializer_list<double> ilist);
    Point& operator= (std::initializer_list<int> ilist);
    ~Point ();

    Point*  SetGlobalIndex (int index);
    int     GetGlobalIndex () const;
    void    ClearListNeighbours ();
    void    AddPointNeighbour (Point* p);
    std::vector <Point *> GetListNeighbours ();
    void    SetListNeighbours (std::vector <Point *>& list);
    void    RemoveThisNeighbourPoint (Point* p);
    void    LinkToCell (Cell* c);
    void    UnlinkToCell (Cell* c);
    std::vector <Cell*> GetLinkedCell () const;
    Point*  DetachFromAll ();

    Point& operator+= (Point&& p);
    Point& operator-= (Point&& p);
    Point& operator*= (Point&& p);
    Point& operator+= (Point& p);
    Point& operator-= (Point& p);
    Point& operator*= (Point& p);
    Point& operator+= (double value);
    Point& operator-= (double value);
    Point& operator*= (double value);
    Point& operator/= (double value);

    std::vector<double> data ();

    double  EuclidianNorm ();
    Point&  Normalize ();

    friend std::ostream & operator<< (std::ostream &out, const Point &p);

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
