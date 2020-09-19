#include <Soliton>
#include "functions.h"

int main(int argc, char *argv[])
{
    HEADERFUN("main");

    Eigen::setNbThreads(4);
    Eigen::initParallel();

    (void)argc;
    (void)argv;

    SolitonRoutines::PrintExecutableName ("POISSON EQUATION - EXECUTABLE CONVERGENCE");
    SolitonRoutines::PrintMacros ();

    InputDatStruct inputdatfile;

    ParseInputDatFile (&inputdatfile, "../input.dat");

    std::vector<double> err_l2, err_l1, err_linf, err_rela, hpres;

    for (double hsz = 0.1; hsz > 3.9*1e-3; hsz-=1e-3)
    {
        hpres.push_back (hsz);
        inputdatfile.hsize = hsz;
        inputdatfile.xm = 0.0;
        inputdatfile.xp = 1.0;
        inputdatfile.ym = 0.0;
        inputdatfile.yp = 1.0;

        std::cout.setstate(std::ios_base::failbit);

        Sto4Sol store;
        SolitonRoutines::Generation (&store, &inputdatfile);


        // ********************************************************************* //
        BEGIN << COLOR_YELLOW << "Solver" << ENDLINE;

        UpgradableSolver<SCHEME_TEMPORAL_SOLVER::NO_TIME> solver(&store);

        solver << ItemSolver<ITEM_SOLVER_TYPE::GRAD_GRAD> ().SetTagToApply (TAG_PHYSICAL::TAG_DOMAIN);
        solver << ItemSolver<ITEM_SOLVER_TYPE::SECOND_MEMBER> ().SetTagToApply (TAG_PHYSICAL::TAG_DOMAIN).SetFunction (fun_secondmember);
        solver << ItemSolver<ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY> ().SetTagToApply (TAG_PHYSICAL::TAG_WALL) .SetFunction (fun_dirichlet);
        solver << ItemSolver<ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY> ().SetTagToApply (TAG_PHYSICAL::TAG_INLET) .SetFunction (fun_dirichlet);
        solver << ItemSolver<ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY> ().SetTagToApply (TAG_PHYSICAL::TAG_OUTLET) .SetFunction (fun_dirichlet);

        solver.SetCoeffPen (1E3);
        solver.SetTime (0.);
        solver.SetTimeStep (1E-4);
        solver.Configure (true);

        PlainVector* sol_num = solver.Compute ();

        std::vector<double> solstd = PlainVector2Vector (sol_num);

        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_num", solstd);

        PlainVector sol_ana;
        FunToVec (&sol_ana, store.mesh, fun_analytic);
        solstd = PlainVector2Vector (&sol_ana);

        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_ana", solstd);

        PlainVector erroabs = GetErrorAbs (store.mesh, &sol_ana, sol_num);
        solstd = PlainVector2Vector (&erroabs);
        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_abs", solstd);

        err_l1.push_back (GetErrorl1 (store.mesh, &sol_ana, sol_num, hsz));
        err_l2.push_back (GetErrorl2 (store.mesh, &sol_ana, sol_num, hsz));
        err_linf.push_back (GetErrorlinf (store.mesh, &sol_ana, sol_num, hsz));
        err_rela.push_back (GetErrorRela (store.mesh, &sol_ana, sol_num, hsz));


        ENDFUN;
        // ********************************************************************* //

        std::cout.clear();


        //        INFOS << COLOR_YELLOW;
        printf ("\033[1;33mdone h = \033[1;32m%10.4e\033[1;33m -- l1-error : \033[1;32m%10.4e\033[1;33m -- l2-error : \033[1;32m%10.4e\033[1;33m -- linf-error : \033[1;32m%10.4e\033[1;33m -- rela-error : \033[1;32m%10.4e\033[1;33m\033[0m\n", hsz, err_l1.back (), err_l2.back (), err_linf.back (), err_rela.back ());
        //        std::cout << ENDLINE;


    }

    std::cout << "err_l1 = [ ";
    for (std::size_t i = 0; i < err_l1.size () - 1; ++i)
        std::cout << err_l1.at (i) << ", ";
    std::cout << err_l1.back () << " ];" << std::endl;

    std::cout << "err_l2 = [ ";
    for (std::size_t i = 0; i < err_l2.size () - 1; ++i)
        std::cout << err_l2.at (i) << ", ";
    std::cout << err_l2.back () << " ];" << std::endl;

    std::cout << "err_linf = [ ";
    for (std::size_t i = 0; i < err_linf.size () - 1; ++i)
        std::cout << err_linf.at (i) << ", ";
    std::cout << err_linf.back () << " ];" << std::endl;

    std::cout << "err_rela = [ ";
    for (std::size_t i = 0; i < err_rela.size () - 1; ++i)
        std::cout << err_rela.at (i) << ", ";
    std::cout << err_rela.back () << " ];" << std::endl;

    std::cout << "h = [ ";
    for (std::size_t i = 0; i < hpres.size () - 1; ++i)
        std::cout << hpres.at (i) << ", ";
    std::cout << hpres.back () << " ];" << std::endl;

    return EXIT_SUCCESS;
}
