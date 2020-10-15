#include <Soliton>
#include "functions.h"

#define NEW_ITEM_WITH_EMB(X)                                                    \
    ItemSolver<ITEM_SOLVER_TYPE::X> ().SetTagToApply (static_cast<PHYS>(item.tagToApply))          \
    .SetDerT (std::size_t(item.dert)).SetVariableOverTime (item.varTime)        \
    .SetEmbMesh (store.listobjects.at (std::size_t(item.id_obj))).SetFunction (fun)

#define NEW_ITEM_WITHOUT_EMB(X)                                                 \
    ItemSolver<ITEM_SOLVER_TYPE::X> ().SetTagToApply (static_cast<PHYS>(item.tagToApply))          \
    .SetDerT (std::size_t(item.dert)).SetVariableOverTime (item.varTime).SetFunction (fun)

#define NEW_ITEM(X)\
    (item.id_obj >= 0 ? NEW_ITEM_WITH_EMB(X) : NEW_ITEM_WITHOUT_EMB(X))

int main(int argc, char *argv[])
{
    HEADERFUN("main");

    Eigen::setNbThreads(4);
    Eigen::initParallel();

    VOID_USE(argc);

    if (argc < 2)
    {
        ERROR << "Use : executable input.slt" << ENDLINE;
        return EXIT_FAILURE;
    }


    VOID_USE(argv);

    SolitonRoutines::PrintExecutableName ("GENERIC EXECUTABLE FOR SOLITON");
    SolitonRoutines::PrintMacros ();

    InputDatStruct inputdatfile;

    ParseInputDatFile (&inputdatfile, argv [1]);

    std::vector<double> err_l2, err_l1, err_linf, err_rela, hpres;

    COUT << std::setw(13) << "hsize" << std::setw(23) << "l1-error" << std::setw(17) << "l2-error" << std::setw(25)<< "linf-error" << std::setw(25) << "rela-error (%)" << ENDLINE;

    for (double hsz = 0.1; hsz > 3.9*1e-3; hsz-=1e-3)
    {
        hpres.push_back (hsz);
        inputdatfile.hsize = hsz;
        inputdatfile.grid_x_m = 0.0;
        inputdatfile.grid_x_p = 1.0;
        inputdatfile.grid_y_m = 0.0;
        inputdatfile.grid_y_p = 1.0;

        std::cout.setstate(std::ios_base::failbit);

        Sto4Sol store;
        SolitonRoutines::Generation (&store, &inputdatfile);


        // ********************************************************************* //
        BEGIN << COLOR_YELLOW << "Solver" << ENDLINE;

        UpgradableSolverBase* solver = NewUpgradableSolver (&store, inputdatfile.scheme);

        for (std::size_t id = 0; id < inputdatfile.items.size (); ++id)
        {
            auto item = inputdatfile.items.at (id);

            std::function<double (Point, double)> fun;

            if (item.fun == "fun_dirichlet")
                fun = fun_dirichlet;
            else if (item.fun == "fun_secondmember")
                fun = fun_secondmember;
            else
                fun = fun_analytic;

            switch (item.type)
            {
            case ITEM_SOLVER_TYPE::EMPTY:
                *solver <<  NEW_ITEM(EMPTY);
                continue;
            case ITEM_SOLVER_TYPE::PHI_PHI:
                *solver <<  NEW_ITEM(PHI_PHI);
                continue;
            case ITEM_SOLVER_TYPE::DAMPING_PHI_PHI:
                *solver <<  NEW_ITEM(DAMPING_PHI_PHI);
                continue;
            case ITEM_SOLVER_TYPE::GRAD_GRAD:
                *solver <<  NEW_ITEM(GRAD_GRAD);
                continue;
            case ITEM_SOLVER_TYPE::SECOND_MEMBER:
                *solver <<  NEW_ITEM(SECOND_MEMBER);
                continue;
            case ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY:
                *solver <<  NEW_ITEM(DIRICHLET_BOUNDARY);
                continue;
            case ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY:
                *solver <<  NEW_ITEM(NEUMANN_BOUNDARY);
                continue;
            case ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM:
                *solver <<  NEW_ITEM(DIRICHLET_BOUNDARY_SBM);
                continue;
            case ITEM_SOLVER_TYPE::NEUMANN_BOUNDARY_SBM:
                *solver <<  NEW_ITEM(NEUMANN_BOUNDARY_SBM);
                continue;
            }
        }

        solver->SetCoeffPen (inputdatfile.penal);
        solver->SetTime (0.);
        solver->SetTimeStep (inputdatfile.dt);
        solver->Configure ();

        PlainVector* sol_num = solver->Compute ();

        std::vector<double> solstd = PlainVector2Vector (sol_num);

        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_num", solstd);

        PlainVector sol_ana;
        FunToVec (&sol_ana, store.mesh, fun_analytic);
        solstd = PlainVector2Vector (&sol_ana);

        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_ana", solstd);

        PlainVector erroabs = GetErrorAbs (store.mesh, &sol_ana, sol_num);
        solstd = PlainVector2Vector (&erroabs);
        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_abs", solstd);

        err_l1.push_back (GetErrorl1 (store.mesh, &sol_ana, sol_num));
        err_l2.push_back (GetErrorl2 (store.mesh, &sol_ana, sol_num));
        err_linf.push_back (GetErrorlinf (store.mesh, &sol_ana, sol_num));
        err_rela.push_back (GetErrorRela (store.mesh, &sol_ana, sol_num));


        ENDFUN;
        // ********************************************************************* //

        std::cout.clear();

        COUT << std::scientific << FLUSHLINE;
        COUT << std::setw(15) << COLOR_YELLOW << hsz << COLOR_DEFAULT << std::setw(15) << COLOR_RED << err_l1.back () << COLOR_DEFAULT << std::setw(12);
        COUT << COLOR_RED << err_l2.back () << COLOR_DEFAULT << std::setw(18) << COLOR_RED << err_linf.back () << COLOR_DEFAULT;
        COUT << std::setw(16) <<  COLOR_RED << err_rela.back () << COLOR_DEFAULT << " (" << COLOR_RED << static_cast<int>(err_rela.back () * 100.) <<  COLOR_DEFAULT << "%)" << ENDLINE;
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
