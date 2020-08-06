#ifndef MATH_H
#define MATH_H

#include <vector>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Core/Defs4Soliton/defs4soliton.h>

class Mesh;
class Point;
class Cell;

namespace std {
template <typename T>
using matrix = std::vector <std::vector <T>>;
}

typedef Eigen::SparseMatrix<double, Eigen::RowMajor>    SparseMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1>        PlainVector;

namespace Math {


double Det (std::matrix<double>& mat);
std::matrix<double> Identity (std::size_t n);
std::matrix<double> Ones (std::size_t n);
std::matrix<double> Zeros (std::size_t n);

}

void FunToVec (PlainVector* out, Mesh * mesh, double (*f) (Point, double), double t = 0.);

void FunToVec (PlainVector* out, Mesh * mesh, double value = 0.);

#endif // MATH_H
