#include "method.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/Eigen>

SOLITON_RETURN Solver::AutoDeduceBest (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display, int maxiter, double eps)
{
    if (A->isApprox (A->adjoint()))
        return Solver::CG (A, b, sol, display, maxiter, eps);
    return Solver::BiCGStab (A, b, sol, display, maxiter, eps);
}

SOLITON_RETURN Solver::CG (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display, int maxiter, double eps)
{
    HEADERFUN("Conjuguate Gradient");

    if (display)
    {
#ifdef VERBOSE
        BEGIN << "conjuguate gradient for symetric matrix." << ENDLINE;
#endif
    }

    //    Eigen::Index n = A->rows ();
    Eigen::Index m = A->cols ();

    if (b->size() != m)
    {
        ERROR << "solver matrix and vector are not the same size." << BLINKRETURN << ENDLINE;
        return SOLITON_FAILURE;
    }

    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower|Eigen::Upper> eigsolver;

    if (maxiter > 0)
        eigsolver.setMaxIterations (maxiter);
    if (eps > 0)
        eigsolver.setTolerance(eps);

    eigsolver.compute(*A);
    *sol = eigsolver.solve(*b);

    if (display)
    {
#ifdef VERBOSE
        INFOS << "iterations        : " << eigsolver.iterations() << ENDLINE;
        INFOS << "estimated error   : " << eigsolver.error() << ENDLINE;
#endif
    }

    return SOLITON_SUCCESS;
}



SOLITON_RETURN Solver::BiCGStab (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display, int maxiter, double eps)
{
    HEADERFUN("Conjuguate Gradient");

    if (display)
    {
#ifdef VERBOSE
        BEGIN << "conjuguate gradient for symetric matrix." << ENDLINE;
#endif
    }

    //    Eigen::Index n = A->rows ();
    Eigen::Index m = A->cols ();

    if (b->size() != m)
    {
        ERROR << "solver matrix and vector are not the same size." << BLINKRETURN << ENDLINE;
        return SOLITON_FAILURE;
    }

    Eigen::BiCGSTAB<SparseMatrix> eigsolver;

    if (maxiter > 0)
        eigsolver.setMaxIterations (maxiter);
    if (eps > 0)
        eigsolver.setTolerance(eps);

    eigsolver.compute(*A);
    *sol = eigsolver.solve(*b);

    if (display)
    {
#ifdef VERBOSE
        INFOS << "iterations        : " << eigsolver.iterations() << ENDLINE;
        INFOS << "estimated error   : " << eigsolver.error() << ENDLINE;
#endif
    }

    return SOLITON_SUCCESS;
}
