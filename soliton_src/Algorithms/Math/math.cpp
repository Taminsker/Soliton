#include "math.h"

#include <Core/core.h>

Matrix3x3 RotateZMatrix (double angle)
{

    Matrix3x3 out;
    out.setZero ();
    out.row (0) << std::cos (angle), -std::sin (angle), 0.;
    out.row (1) << std::sin (angle), std::cos (angle), 0.;
    out.row (2) << 0., 0., 1.;
    return out;
}
//double Det (unstd::matrix<double>& mat)
//{

//    std::size_t n = mat.size ();
//    double det = 0;

//    if (n != mat[0].size ())
//    {
//        ERROR << "call on a matrix which is not square" << BLINKRETURN << ENDLINE;
//        return 0.0;
//    }

//    if (n == 1)
//        return mat[0][0];

//    for (std::size_t id = 0; id < n; ++id)
//    {
//        unstd::matrix<double> matcopy = mat;

//        for (std::size_t i = 0; i < n; ++i)
//            for (std::size_t j = 0; j < n; ++j)
//                if (j == id)
//                    matcopy[i].erase (matcopy[i].begin () + static_cast<int>(id));

//        matcopy.erase (matcopy.begin ());

//        det += std::pow (-1, static_cast<int>(id)) * mat[0][id] * Math::Det (matcopy);
//    }

//    return det;
//}

//unstd::matrix<double> Math::Identity (std::size_t n)
//{
//    unstd::matrix<double> mat(n, std::vector<double>(n));
//    for (std::size_t id = 0; id < n; id++)
//        mat [id][id] = 1.;

//    return mat;
//}

//unstd::matrix<double> Math::Ones (std::size_t n)
//{
//    return unstd::matrix<double>(n, std::vector<double>(n, 1.));
//}

//unstd::matrix<double> Math::Zeros (std::size_t n)
//{
//    return unstd::matrix<double> (n, std::vector<double>(n, 0.));
//}

//unstd::matrix<double> Math::Transpose (const unstd::matrix<double>& mat)
//{
//    unstd::matrix<double> out = mat;
//    Math::TransposeInPlace (&out);

//    return out;
//}


//void Math::TransposeInPlace (unstd::matrix<double>* mat)
//{
//    std::size_t n = mat->size ();

//    for (std::size_t i = 0; i < n; ++i)
//        for (std::size_t j = 0; j < i; ++j)
//        {
//            double value = mat->operator[] (i)[j];
//            mat->at (i).at (j) = mat->at (j).at (i);
//            mat->at (j).at (i) = value;
//        }

//    return;
//}

std::vector<double> PlainVector2Vector (PlainVector* vec)
{
    return std::vector<double>(vec->data(), vec->data() + vec->rows() * vec->cols());
}


void FunToVec (PlainVector* out, Mesh * mesh, std::function<double(Point, double)> f, double t)
{
    int n = mesh->GetNumberOfPoints ();

    out->resize (n);
    out->setZero ();

    for (int i = 0; i < n; i++)
        out->coeffRef (i) = f(*mesh->GetPoint (i), t);

    return;
}

void FunToVec (PlainVector* out, Mesh * mesh, double value)
{
    out->setConstant (mesh->GetNumberOfPoints (), value);
    return;
}
