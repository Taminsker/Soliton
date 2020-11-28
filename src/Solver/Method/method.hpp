#ifndef METHOD_H
#define METHOD_H

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Core/core.h>
#include <Algorithms/Math/math.h>

namespace Solver {
void AutoDeduceBest (const SparseMatrix_eig *A, const PlainVector_eig *b, PlainVector_eig *sol, bool display = true, int maxiter = -1, real_t eps = -1.);

void CG (const SparseMatrix_eig *A, const PlainVector_eig *b, PlainVector_eig *sol, bool display = true, int maxiter = -1, real_t eps = -1.);

void BiCGStab (const SparseMatrix_eig *A, const PlainVector_eig *b, PlainVector_eig *sol, bool display = true, int maxiter = -1, real_t eps = -1.);

}

void ImposeDirichlet (Mesh* mesh, SparseMatrix_eig *A, PlainVector_eig *secondMember,
                      real_t (*g) (Point, real_t), std::vector <int>* listIndex,
                      real_t t = 0.);

real_t          GetErrorl1 (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view = true);
real_t          GetErrorl2 (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view = true);
real_t          GetErrorlinf (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view = true);
PlainVector_eig GetErrorAbs (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view = true);
PlainVector_eig GetErrorRelaPercent (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view = true);
real_t          GetErrorRela (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view = true);

#endif // METHOD_H
