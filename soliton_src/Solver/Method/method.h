#ifndef METHOD_H
#define METHOD_H

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Core/core.h>
#include <Algorithms/Math/math.h>

namespace Solver {
SOLITON_RETURN AutoDeduceBest (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display = true, int maxiter = -1, double eps = -1.);

SOLITON_RETURN CG (const SparseMatrix *A, const PlainVector *b, PlainVector* sol,  bool display = true, int maxiter = -1, double eps = -1.);

SOLITON_RETURN BiCGStab (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display = true, int maxiter = -1, double eps = -1.);

}





#endif // METHOD_H
