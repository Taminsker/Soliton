#ifndef METHOD_H
#define METHOD_H

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Core/core.h>
#include <Algorithms/Math/math.h>

namespace Solver {
void AutoDeduceBest (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display = true, int maxiter = -1, double eps = -1.);

void CG (const SparseMatrix *A, const PlainVector *b, PlainVector* sol,  bool display = true, int maxiter = -1, double eps = -1.);

void BiCGStab (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display = true, int maxiter = -1, double eps = -1.);

}

void ImposeDirichlet (Mesh* mesh, SparseMatrix* A, PlainVector* secondMember,
                      double (*g) (Point, double), std::vector <int>* listIndex,
                      double t = 0.);

double      GetErrorl1 (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num);
double      GetErrorl2 (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num);
double      GetErrorlinf (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num);
PlainVector GetErrorAbs (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num);
PlainVector GetErrorRelaPercent (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num);
double      GetErrorRela (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num);

#endif // METHOD_H
