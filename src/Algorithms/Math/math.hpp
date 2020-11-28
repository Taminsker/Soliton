#ifndef MATH_H
#define MATH_H

#include <vector>
#include <cmath>

#include <ProgDef/proddef.h>

#include "Iterator/iterators.h"

class Mesh;
class Point;
class Cell;
class SolitonFunctor;

class Duo
{
public:
    Duo () : m_id (0), m_value(0) {}

    Duo (const int& id, const real_t& v = real_t(0))
        : m_id(id), m_value(v)
    {}

    int idx () const { return m_id;}
    const real_t& value() const { return m_value; }

protected:
    int m_id;
    real_t m_value;
};

Matrix3x3_eig    RotateZMatrix (real_t angle);
std::vector<real_t> PlainVector2Vector (PlainVector_eig *vec);
void        FunToVec (PlainVector_eig *out, Mesh * mesh, std::function<real_t(Point, real_t)> f, real_t t = 0.);
void        FunToVec (PlainVector_eig *out, Mesh * mesh, SolitonFunctor* f, real_t t = 0.);
void        FunToVec (PlainVector_eig *out, Mesh * mesh, real_t value = 0.);

#endif // MATH_H
