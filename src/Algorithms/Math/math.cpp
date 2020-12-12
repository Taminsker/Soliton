#include "math.h"

#include "../../Core/core.hpp"
#include "../../Solver/SolitonFunctor/solitonfunctor.hpp"

Matrix3x3
RotateZMatrix (real_t angle)
{
    Matrix3x3 out;
    out.setZero ();
    out.row (0) << std::cos (angle), -std::sin (angle), 0.;
    out.row (1) << std::sin (angle), std::cos (angle), 0.;
    out.row (2) << 0., 0., 1.;
    return out;
}
//real_t Det (unstd::matrix<real_t>& mat)
//{

//  ul_t n = mat.size ();
//  real_t det = 0;

//  if (n != mat[0].size ())
//  {
//    ERROR << "call on a matrix which is not square" << BLINKRETURN << ENDLINE;
//    return 0.0;
//  }

//  if (n == 1)
//    return mat[0][0];

//  for (ul_t id = 0; id < n; ++id)
//  {
//    unstd::matrix<real_t> matcopy = mat;

//    for (ul_t i = 0; i < n; ++i)
//      for (ul_t j = 0; j < n; ++j)
//        if (j == id)
//          matcopy[i].erase (matcopy[i].begin () + static_cast<int>(id));

//    matcopy.erase (matcopy.begin ());

//    det += std::pow (-1, static_cast<int>(id)) * mat[0][id] * Math::Det (matcopy);
//  }

//  return det;
//}

//unstd::matrix<real_t> Math::Identity (ul_t n)
//{
//  unstd::matrix<real_t> mat(n, std::vector<real_t>(n));
//  for (ul_t id = 0; id < n; id++)
//    mat [id][id] = 1.;

//  return mat;
//}

//unstd::matrix<real_t> Math::Ones (ul_t n)
//{
//  return unstd::matrix<real_t>(n, std::vector<real_t>(n, 1.));
//}

//unstd::matrix<real_t> Math::Zeros (ul_t n)
//{
//  return unstd::matrix<real_t> (n, std::vector<real_t>(n, 0.));
//}

//unstd::matrix<real_t> Math::Transpose (const unstd::matrix<real_t>& mat)
//{
//  unstd::matrix<real_t> out = mat;
//  Math::TransposeInPlace (&out);

//  return out;
//}

//void Math::TransposeInPlace (unstd::matrix<real_t>* mat)
//{
//  ul_t n = mat->size ();

//  for (ul_t i = 0; i < n; ++i)
//    for (ul_t j = 0; j < i; ++j)
//    {
//      real_t value = mat->operator[] (i)[j];
//      mat->at (i).at (j) = mat->at (j).at (i);
//      mat->at (j).at (i) = value;
//    }

//  return;
//}

std::vector<real_t>
PlainVector2Vector (DenseVector * vec)
{
    return std::vector<real_t> (vec->data (), vec->data () + vec->rows () * vec->cols ());
}

void
FunToVec (DenseVector * out, Mesh * mesh, std::function<real_t (Point, real_t)> f, real_t t)
{
    int n = mesh->GetNumberOfPoints ();

    out->resize (n);
    out->setZero ();

    for (int i = 0; i < n; i++)
        out->coeffRef (i) = f (*mesh->GetPoint (i), t);

    return;
}

void
FunToVec (DenseVector * out, Mesh * mesh, SolitonFunctor * f, real_t t)
{
    int n = mesh->GetNumberOfPoints ();

    out->resize (n);
    out->setZero ();

    for (int i = 0; i < n; i++)
        out->coeffRef (i) = f->ToReal (mesh->GetPoint (i), t, nullptr);

    return;
}

void
FunToVec (DenseVector * out, Mesh * mesh, real_t value)
{
    out->setConstant (mesh->GetNumberOfPoints (), value);
    return;
}
