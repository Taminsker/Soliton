#ifndef MATH_H
#define MATH_H

#include <vector>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <ProgDef/proddef.h>

#include "Iterator/iterators.h"

class Mesh;
class Point;
class Cell;

namespace unstd {
template <typename T>
using matrix = std::vector <std::vector <T>>;
}

typedef Eigen::SparseMatrix<double, Eigen::RowMajor>    SparseMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1>        PlainVector;
typedef Eigen::Matrix<double, 3, 3>                     Matrix3x3;
typedef Eigen::Triplet<double>                          Triplet;

class Duo
{
public:
  Duo () : m_id (0), m_value(0) {}

  Duo (const int& id, const double& v = double(0))
    : m_id(id), m_value(v)
  {}

  int idx () const { return m_id;}
  const double& value() const { return m_value; }

protected:
  int m_id;
  double m_value;
};

Matrix3x3 RotateZMatrix (double angle);
std::vector<double> PlainVector2Vector (PlainVector* vec);
void FunToVec (PlainVector* out, Mesh * mesh, std::function<double(Point, double)> f, double t = 0.);
void FunToVec (PlainVector* out, Mesh * mesh, double value = 0.);

#endif // MATH_H
