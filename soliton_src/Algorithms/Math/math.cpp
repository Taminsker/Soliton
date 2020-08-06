#include "math.h"

#include <Core/core.h>

double Math::Det (std::matrix<double>& mat)
{
    HEADERFUN("ComputeDet");

    std::size_t n = mat.size ();
    double det = 0;

    if (n != mat[0].size ())
    {
        ERROR << "call on a matrix which is not square" << BLINKRETURN << ENDLINE;
        return 0.0;
    }

    if (n == 1)
        return mat[0][0];

    for (std::size_t id = 0; id < n; ++id)
    {
        std::matrix<double> matcopy = mat;

        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                if (j == id)
                    matcopy[i].erase (matcopy[i].begin () + int(id));

        matcopy.erase (matcopy.begin ());

        det += std::pow (-1, int(id)) * mat[0][id] * Math::Det (matcopy);
    }

    return det;
}

std::matrix<double> Math::Identity (std::size_t n)
{
    std::matrix<double> mat(n, std::vector<double>(n));
    for (std::size_t id = 0; id < n; id++)
        mat [id][id] = 1.;

    return mat;
}

std::matrix<double> Math::Ones (std::size_t n)
{
    return std::matrix<double>(n, std::vector<double>(n, 1.));
}

std::matrix<double> Math::Zeros (std::size_t n)
{
    return std::matrix<double> (n, std::vector<double>(n, 0.));
}

void FunToVec (PlainVector* out, Mesh * mesh, double (*f) (Point, double), double t)
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
