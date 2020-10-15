#include "method.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/Eigen>

void Solver::AutoDeduceBest (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display, int maxiter, double eps)
{
    if (A->isApprox (A->adjoint()))
        return Solver::CG (A, b, sol, display, maxiter, eps);
    return Solver::BiCGStab (A, b, sol, display, maxiter, eps);
}

void Solver::CG (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display, int maxiter, double eps)
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
        return;
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

    return;
}



void Solver::BiCGStab (const SparseMatrix *A, const PlainVector *b, PlainVector* sol, bool display, int maxiter, double eps)
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
        return;
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

    return;
}

void ImposeDirichlet (Mesh* mesh, SparseMatrix* A, PlainVector* secondMember,
                      double (*g) (Point, double), std::vector <int>* listIndex,
                      double t)
{
    HEADERFUN ("ImposeDirichlet");

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

double GetErrorl1 (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num)
{
    VOID_USE(mesh);

    double error_l1;
    PlainVector u_abs = mesh->GetPrescribedSize () * (*u_ana - *u_num).cwiseAbs();
    error_l1 = u_abs.sum();

    INFOS << COLOR_GREEN << "error l^1 has been calculated : " << std::scientific << error_l1 << ENDLINE;
    return error_l1;
}

double GetErrorl2 (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num)
{
    VOID_USE(mesh);

    double error_l2 = mesh->GetPrescribedSize () *(*u_ana - *u_num).norm ();

    INFOS << COLOR_GREEN << "error l^2 has been calculated : " << std::scientific << error_l2 << ENDLINE;
    return error_l2;
}

double GetErrorlinf (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num)
{
    VOID_USE(mesh);

    PlainVector u_abs = mesh->GetPrescribedSize () * (*u_ana - *u_num).cwiseAbs();
    double error_linf = u_abs.maxCoeff();

    INFOS << COLOR_GREEN << "error l^inf has been calculated : " << std::scientific << error_linf << ENDLINE;

    return error_linf;
}

PlainVector GetErrorAbs (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num)
{
    VOID_USE(mesh);
    INFOS << COLOR_GREEN << "error abs has been calculated." << ENDLINE;

    return (*u_ana - *u_num).cwiseAbs ();
}

PlainVector GetErrorRelaPercent (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num)
{
    VOID_USE(mesh);
    INFOS << COLOR_GREEN << "error rela percent has been calculated." << ENDLINE;

    return  (*u_ana - *u_num).cwiseAbs ().cwiseQuotient((u_ana->cwiseAbs ().array () + EPSILON).matrix ());
}

double GetErrorRela (Mesh * mesh, const PlainVector *u_ana, const PlainVector *u_num)
{
    VOID_USE(mesh);

    double error_rela = (*u_ana - *u_num).norm () / u_ana->norm ();
    INFOS << COLOR_GREEN << "error rela has been calculated : " << std::scientific << error_rela << std::defaultfloat << " " << error_rela * 100. << "%..." << ENDLINE;

    return error_rela;
}
