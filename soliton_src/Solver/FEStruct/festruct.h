#ifndef FESTRUCT_H
#define FESTRUCT_H

#include <vector>
#include <string>
#include <functional>

#include <Core/core.h>
#include <Algorithms/Math/math.h>

#define LAMBDA_EVAL(x)\
    [] (Point p) -> double {(void)p; return x;}

#define LAMBDA_FUNCTION(x)\
    [] (Point p) -> double {(void)p; x}

struct Functor_L2P;

// This part is a free adaptation of the list of elements present in
// 'Méthode des éléments finis' de Gouri Dhatt et al.

typedef enum
{
    LAGRANGE,
    END_NAME_FE
} NAME_FE;

template <NAME_FE namefe, VTKCellType T>
class FEElement
{
public:
    typedef std::vector <std::function <double (Point)>> FuncContainer;

    FEElement () :
        name (T),
        _namefe (namefe),
        npts (0),
        coor ({}),
        fun ({}),
        der ({}),
        volRef (0)
    {
        Build ();
    }

    ~FEElement ()
    {
        name = VTK_EMPTY_CELL;
        npts = 0;
        coor.clear ();
        fun.clear ();
        der.clear ();
    }

    std::matrix<double> Jac (Cell* cell)
    {
        (void)cell;
        return Math::Zeros(3);
    }

    Point B (Cell*cell)
    {
        (void)cell;
        return Point();
    }

    Point Trans (Cell* cell, Point* p)
    {
        Point out;
        std::matrix<double> A = Jac(cell);
        Point B = this->B(cell);

        out.x = A[0][0] * p->x + A[0][1] * p->y + A[0][2] * p->z + B.x;
        out.y = A[1][0] * p->x + A[1][1] * p->y + A[1][2] * p->z + B.y;
        out.z = A[2][0] * p->x + A[2][1] * p->y + A[2][2] * p->z + B.z;

        return out;
    }

    VTKCellType name;
    NAME_FE _namefe;
    size_t npts;
    std::vector <Point> coor;
    FuncContainer fun;
    std::vector <Functor_L2P> der;
    double volRef;

protected:
    void Build () {}
};

template <VTKCellType T>
using LagrangeFE = FEElement<LAGRANGE, T>;

struct Functor_L2P
{
    Point operator() (Point p)
    {
        Point out;
        out.x = x(p);
        out.y = y(p);
        out.z = z(p);

        return out;
    }

    template <NAME_FE namefe, VTKCellType T>
    friend class FEElement;

protected:
    typedef std::function <double (Point)> Fun;

    Fun x = LAMBDA_EVAL(0);
    Fun y = LAMBDA_EVAL(0);
    Fun z = LAMBDA_EVAL(0);
};

#define NAMELAGRANGEFE(x) _##x##_
#define DECLLAGRANGEFE(x) LagrangeFE<x>* NAMELAGRANGEFE(x)

class FEStruct
{

public:
    FEStruct ();
    ~FEStruct ();

    template<VTKCellType T>
    LagrangeFE<T>* GetElementFor ()
    {
        return nullptr;
    }

private:
    DECLLAGRANGEFE (VTK_EMPTY_CELL);
    DECLLAGRANGEFE (VTK_LINE);
    DECLLAGRANGEFE (VTK_VERTEX);
    DECLLAGRANGEFE (VTK_TRIANGLE);
    DECLLAGRANGEFE (VTK_QUADRATIC_EDGE);
    DECLLAGRANGEFE (VTK_QUADRATIC_TRIANGLE);

};

#endif // FESTRUCT_H
