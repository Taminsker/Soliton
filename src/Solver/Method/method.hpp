#ifndef SRC_SOLVER_METHOD_METHOD_HPP
#define SRC_SOLVER_METHOD_METHOD_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../../Algorithms/Math/math.hpp"
#include "../../Core/core.hpp"

namespace Solver
{
void AutoDeduceBest (const SparseMatrix * A, const DenseVector * b, DenseVector * sol, bool display = true, int maxiter = -1, real_t eps = -1.);

void CG (const SparseMatrix * A, const DenseVector * b, DenseVector * sol, bool display = true, int maxiter = -1, real_t eps = -1.);

void BiCGStab (const SparseMatrix * A, const DenseVector * b, DenseVector * sol, bool display = true, int maxiter = -1, real_t eps = -1.);

}  // namespace Solver

void ImposeDirichlet (Mesh * mesh, SparseMatrix * A, DenseVector * secondMember,
                      real_t (*g) (Point, real_t), std::vector<int> * listIndex,
                      real_t t = 0.);

real_t      GetErrorl1 (Mesh * mesh, const DenseVector * u_ana, const DenseVector * u_num, bool view = true);
real_t      GetErrorl2 (Mesh * mesh, const DenseVector * u_ana, const DenseVector * u_num, bool view = true);
real_t      GetErrorlinf (Mesh * mesh, const DenseVector * u_ana, const DenseVector * u_num, bool view = true);
DenseVector GetErrorAbs (Mesh * mesh, const DenseVector * u_ana, const DenseVector * u_num, bool view = true);
DenseVector GetErrorRelaPercent (Mesh * mesh, const DenseVector * u_ana, const DenseVector * u_num, bool view = true);
real_t      GetErrorRela (Mesh * mesh, const DenseVector * u_ana, const DenseVector * u_num, bool view = true);

#endif /* SRC_SOLVER_METHOD_METHOD_HPP */
