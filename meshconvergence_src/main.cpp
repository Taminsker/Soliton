#include <Soliton>
#include "functions.h"

int main(int argc, char *argv[])
{
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

    for (double hsz = 0.1; hsz > 9.9*1e-3; hsz-=1e-3)
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

        SuperSolver solver (&store, inputdatfile.listTagsSolver);

        for (ItemSolverDatStruct itemdat : inputdatfile.items)
        {
            SuperItem item (itemdat.type);
            item.SetTag2Apply (itemdat.tagToApply);
            item.SetDerT(static_cast<std::size_t>(itemdat.dert));

            if (itemdat.id_obj > 0)
                item.SetTargetMesh (&store.listobjects.at (static_cast<std::size_t> (itemdat.id_obj)));

            if (itemdat.fun == "fun_dirichlet")
                item.SetFunction2Apply (fun_dirichlet);
            else if (itemdat.fun == "fun_neumann")
                item.SetFunction2Apply (fun_neumann);
            else if (itemdat.fun == "fun_secondmember")
                item.SetFunction2Apply (fun_secondmember);
            else if (itemdat.fun == "fun_analytic")
                item.SetFunction2Apply (fun_analytic);
            else
                item.SetFunction2Apply (NULL_FUNC_NAME);


            solver.AddItem (item, static_cast<std::size_t> (itemdat.dert), itemdat.solver_tag);
        }

        solver.SetScheme (inputdatfile.scheme);
        solver.SetCoeffPen (inputdatfile.penal);
        solver.SetTime (0.);
        solver.SetTimeStep (inputdatfile.dt);
        solver.SetCollectContribution (inputdatfile.colConcItem);

        solver.SetViewOn ();
        solver.SetForceComputeOn ();
        solver.SetCollectContributionOff ();

        solver.Configure ();

        solver.SetViewOff ();
        solver.SetForceComputeOff ();
        solver.SetCollectContributionOff ();

        std::vector<PlainVector*> sols = solver.Compute ();
        std::vector<INTER> list_tags_inter = solver.GetListOfInterTags ();

        //    std::size_t numPoints = static_cast<std::size_t>(store.mesh->GetNumberOfPoints ());
        //    std::size_t numCells = static_cast<std::size_t>(store.mesh->GetNumberOfCells ());
        std::size_t numEdges = static_cast<std::size_t>(store.mesh->GetNumberOfEdges ());

        HetInt::Array* tagSurEdge = store.mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_INTERSECTION);

        double error_l1 = 0.;
        double error_l2 = 0.;
        double error_rela = 0.;
        double error_linf = 0.;

        for (std::size_t idvec = 0; idvec < sols.size (); ++idvec)
        {
            STATUS << "postprocess on solver with tag " <<  COLOR_RED << to_string(list_tags_inter.at (idvec)) << ENDLINE;

            PlainVector* solnum = sols.at (idvec);

            PlainVector sol_ana;
            FunToVec (&sol_ana, store.mesh, fun_analytic);

            if (list_tags_inter.at (idvec) != INTER::DEFAULT)
                for (std::size_t idEdge = 0; idEdge < numEdges; ++idEdge)
                {
                    if (tagSurEdge->vec.at (idEdge) == int(list_tags_inter.at (idvec)) ||
                            tagSurEdge->vec.at (idEdge) == int(INTER::MIXED))
                        continue;

                    for (Point* pt : *store.mesh->GetEdge (int(idEdge))->GetPoints ())
                    {
                        bool pass = false;
                        for (Cell* cell : pt->GetLinkedCell ())
                        {
                            if (cell->GetCat () == CAT_CELL_EDGE::CELL)
                                continue;

                            if (tagSurEdge->vec.at (static_cast<std::size_t> (cell->GetGlobalIndex ())) == int(INTER::MIXED))
                            {
                                pass = true;
                                break;
                            }
                        }

                        if (pass)
                            continue;

                        solnum->coeffRef (pt->GetGlobalIndex ()) = -10.;
                        sol_ana.coeffRef (pt->GetGlobalIndex ()) = -10.;
                    }
                }


            std::string tagsolver = "_" + to_string(list_tags_inter.at (idvec));

            if (tagsolver ==  "_" + to_string(INTER::DEFAULT))
                tagsolver = "";

            std::vector<double> solstd = PlainVector2Vector (solnum);
            store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_num" + tagsolver, solstd);

            solstd = PlainVector2Vector (&sol_ana);
            store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_ana" + tagsolver, solstd);

            PlainVector erroabs = GetErrorAbs (store.mesh, &sol_ana, solnum);
            solstd = PlainVector2Vector (&erroabs);
            store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_abs" + tagsolver, solstd);

    //        PlainVector errorelapercent = GetErrorRelaPercent (store.mesh, &sol_ana, solnum);
    //        solstd = PlainVector2Vector (&errorelapercent);
    //        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_rela_percent" + tagsolver, solstd);

            error_l1 = GetErrorl1     (store.mesh, &sol_ana, solnum);
            error_l2 = GetErrorl2 (store.mesh, &sol_ana, solnum);
            error_linf = GetErrorlinf (store.mesh, &sol_ana, solnum);
            error_rela = GetErrorRela (store.mesh, &sol_ana, solnum);

        }

        err_l1.push_back    (error_l1);
        err_l2.push_back    (error_l2);
        err_linf.push_back  (error_linf);
        err_rela.push_back  (error_rela);

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
