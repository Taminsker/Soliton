#ifndef SRC_ALGORITHMS_MATH_MATH_HPP
#define SRC_ALGORITHMS_MATH_MATH_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <vector>

#include "../../solitonheader.hpp"
#include "Iterator/iterators.hpp"

// Math def

namespace unstd
{
template <typename T>
using matrix = std::vector<std::vector<T>>;
}

typedef Eigen::SparseMatrix<real_t, Eigen::RowMajor> SparseMatrix;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1>     DenseVector;
typedef Eigen::Matrix<real_t, 3, 3>                  Matrix3x3;
typedef Eigen::Triplet<real_t>                       Triplet;

class Mesh;
class Point;
class Cell;
class SolitonFunctor;

class Duo
{
public:
    Duo () : m_id (0),
             m_value (0) {}

    Duo (const int & id, const real_t & v = real_t (0))
        : m_id (id),
          m_value (v)
    {
    }

    int
    idx () const
    {
        return m_id;
    }
    const real_t &
    value () const
    {
        return m_value;
    }

protected:
    int    m_id;
    real_t m_value;
};

Matrix3x3           RotateZMatrix (real_t angle);
std::vector<real_t> PlainVector2Vector (DenseVector * vec);
void                FunToVec (DenseVector * out, Mesh * mesh, std::function<real_t (Point, real_t)> f, real_t t = 0.);
void                FunToVec (DenseVector * out, Mesh * mesh, SolitonFunctor * f, real_t t = 0.);
void                FunToVec (DenseVector * out, Mesh * mesh, real_t value = 0.);

#endif /* SRC_ALGORITHMS_MATH_MATH_HPP */
