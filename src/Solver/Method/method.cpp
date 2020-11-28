#include "method.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/Eigen>

void
Solver::AutoDeduceBest (const SparseMatrix_eig *A, const PlainVector_eig *b,
                        PlainVector_eig* sol, bool display, int maxiter, real_t eps)
{

    if (A->nonZeros () == 0 || *b == PlainVector_eig::Zero (b->size ()))
    {
        sol->setZero (b->size ());
        return;
    }

    if (A->isApprox (A->adjoint()))
        return Solver::CG (A, b, sol, display, maxiter, eps);
    return Solver::BiCGStab (A, b, sol, display, maxiter, eps);
}

void
Solver::CG (const SparseMatrix_eig *A, const PlainVector_eig *b,
            PlainVector_eig* sol, bool display, int maxiter, real_t eps)
{
    if (display)
    {
        //#ifdef VERBOSE
        //    BEGIN << "conjuguate gradient for symetric matrix." << ENDLINE;
        //#endif
    }

    //  Eigen::Index n = A->rows ();
    Eigen::Index m = A->cols ();

    if (b->size() != m)
    {
        ERROR << "solver matrix and vector are not the same size." << BLINKRETURN << ENDLINE;
        return;
    }

    Eigen::ConjugateGradient<SparseMatrix_eig, Eigen::Lower|Eigen::Upper> eigsolver;

    if (maxiter > 0)
        eigsolver.setMaxIterations (maxiter);
    if (eps > 0)
        eigsolver.setTolerance(eps);

    eigsolver.compute(*A);
    *sol = eigsolver.solve(*b);

    if (display)
    {
#ifdef VERBOSE
        INFOS << "iterations    : " << eigsolver.iterations() << ENDLINE;
        INFOS << "estimated error  : " << eigsolver.error() << ENDLINE;
        ENDFUN;
#endif
    }

    return;
}



void
Solver::BiCGStab (const SparseMatrix_eig *A, const PlainVector_eig *b,
                  PlainVector_eig* sol, bool display, int maxiter, real_t eps)
{
    if (display)
    {
        //#ifdef VERBOSE
        //    BEGIN << "conjuguate gradient for symetric matrix." << ENDLINE;
        //#endif
    }

    //  Eigen::Index n = A->rows ();
    Eigen::Index m = A->cols ();

    if (b->size() != m)
    {
        ERROR << "solver matrix and vector are not the same size." << BLINKRETURN << ENDLINE;
        return;
    }

    Eigen::BiCGSTAB<SparseMatrix_eig> eigsolver;

    if (maxiter > 0)
        eigsolver.setMaxIterations (maxiter);
    if (eps > 0)
        eigsolver.setTolerance(eps);

    eigsolver.compute(*A);
    *sol = eigsolver.solve(*b);

    if (display)
    {
#ifdef VERBOSE
        INFOS << "iterations    : " << eigsolver.iterations() << ENDLINE;
        INFOS << "estimated error  : " << eigsolver.error() << ENDLINE;
        ENDFUN;
#endif
    }


    return;
}

void
ImposeDirichlet (Mesh* mesh, SparseMatrix_eig* A, PlainVector_eig* secondMember,
                 real_t (*g) (Point, real_t), std::vector <int>* listIndex,
                 real_t t)
{
    INFOS << "impose Dirichlet on " << listIndex->size () << " points." << ENDLINE;

    for (int i : *listIndex)
        A->row (i) *= 0.;

    *A = A->transpose (); // Pour pouvoir manipuler les colonnes de A

    for (int i : *listIndex) {

        *secondMember -= g(*mesh->GetPoint (i), t) * A->row (i).transpose ();
        A->row (i) *= 0.;
        A->coeffRef (i,i) = 1.;
        secondMember->coeffRef (i) = g(*mesh->GetPoint (i), t);
    }

    *A = A->transpose ().pruned ();

    return;
}

real_t
GetErrorl1 (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view)
{
    real_t error_l1;
    PlainVector_eig u_abs = mesh->GetPrescribedSize () * (*u_ana - *u_num).cwiseAbs();
    error_l1 = u_abs.sum();

    if (view)
        INFOS << COLOR_GREEN << "error l^1   : " << std::scientific << error_l1 << ENDLINE;
    return error_l1;
}

real_t
GetErrorl2 (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view)
{
    real_t error_l2 = mesh->GetPrescribedSize () *(*u_ana - *u_num).norm ();

    if (view)
        INFOS << COLOR_GREEN << "error l^2   : " << std::scientific << error_l2 << ENDLINE;
    return error_l2;
}

real_t
GetErrorlinf (Mesh * mesh, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view)
{
    PlainVector_eig u_abs = mesh->GetPrescribedSize () * (*u_ana - *u_num).cwiseAbs();
    real_t error_linf = u_abs.maxCoeff();

    if (view)
        INFOS << COLOR_GREEN << "error l^inf : " << std::scientific << error_linf << ENDLINE;

    return error_linf;
}

PlainVector_eig
GetErrorAbs (Mesh *, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view)
{
    if (view)
        INFOS << COLOR_GREEN << "error abs." << ENDLINE;

    return (*u_ana - *u_num).cwiseAbs ();
}

PlainVector_eig
GetErrorRelaPercent (Mesh *, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view)
{
    if (view)
        INFOS << COLOR_GREEN << "error rela percent." << ENDLINE;

    return (*u_ana - *u_num).cwiseAbs ().cwiseQuotient((u_ana->cwiseAbs ().array () + EPSILON).matrix ());
}

real_t
GetErrorRela (Mesh *, const PlainVector_eig *u_ana, const PlainVector_eig *u_num, bool view)
{
    real_t error_rela = (*u_ana - *u_num).norm () / u_ana->norm ();

    if (view)
        INFOS << COLOR_GREEN << "error rela  : " << std::scientific << error_rela << std::defaultfloat << "\t" << COLOR_RED << error_rela * 100. << "%..." << ENDLINE;

    return error_rela;
}