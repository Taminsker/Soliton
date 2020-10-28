#include <Soliton>

double fun_analytic (Point p, double t)
{
    VOID_USE(t);
    return (p.x - p.x * p.x) * (p.y - p.y * p.y);
    //    return 1. * p.x * p.x + 1. * p.y * p.y + 1. * p.x * p.y + 1. * p.x + 1. * p.y + 1.;
    //        return p.x + p.y;
    //    return -0.25 * p.x * p.x - 0.25 * p.y * p.y;
    //    return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
    //        return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
}

double fun_neumann (Point p, double t)
{
    VOID_USE(p);
    VOID_USE(t);

    return 0.;
    //        return (p.x * p.x + p.y * p.y);
}

double fun_dirichlet (Point p, double t)
{
    VOID_USE(p);
    VOID_USE(t);

    //    INFOS << "call dirichlet" << ENDLINE;
    return fun_analytic (p, t);
    //        return 0.;
}

double fun_secondmember (Point p, double t)
{
    VOID_USE(t);
    VOID_USE(p);

    //    return 1.;
    //        return 0;
    return 2. * (p.x - p.x * p.x + p.y - p.y * p.y);
    //    return -2 * (1. + 1.);
    //        return 2. * M_PI * M_PI * fun_analytic (p, t);
}


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

    Sto4Sol store;
    InputDatStruct inputdatfile;

    ParseInputDatFile (&inputdatfile, argv [1]);

    SolitonRoutines::Generation (&store, &inputdatfile);
    SuperSolver solver (&store, inputdatfile.listTagsSolver);

    for (ItemSolverDatStruct itemdat : inputdatfile.items)
    {
        SuperItem item (itemdat.type);
        item.SetTag2Apply (itemdat.tagToApply);
        item.SetDerT(static_cast<std::size_t>(itemdat.dert));
        item.SetVariableOverTime (itemdat.varTime);

        if (itemdat.id_obj >= 0)
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
    solver.SetPowPen (inputdatfile.powpenalty);
    solver.SetTime (0.);
    solver.SetTimeStep (inputdatfile.dt);


    solver.SetCollectContribution (inputdatfile.colConcItem);

    solver.SetViewOn ();
    solver.SetForceComputeOn ();

    solver.Configure ();

    solver.SetViewOff ();
    solver.SetForceComputeOff ();

    std::vector<PlainVector*> list_sols_num = solver.Compute ();

    SolitonRoutines::PostProcess (&store, &solver, &list_sols_num, fun_analytic);

    ENDFUN;

    // ********************************************************************* //
    //
    //    std::vector<std::vector<Point*>> ref(3);
    //
    //    for (int i = 0; i < store.mesh->GetNumberOfCells (); ++i)
    //    {
    //        Point* centroid = store.mesh->GetCell (i)->GetCentroid ();
    //
    //        std::size_t j = 0;
    //        for (Edge *edge : *store.mesh->GetCell (i)->GetEdges ())
    //        {
    //            Point* cent = edge->GetCentroid ();
    //            ref [j].push_back (new Point(*cent - *centroid));
    //            j++;
    //        }
    //    }
    //
    //    store.mesh->GetCellsData ()->GetVecArrays ()->Add ("ref0", ref[0]);
    //    store.mesh->GetCellsData ()->GetVecArrays ()->Add ("ref1", ref[1]);
    //    store.mesh->GetCellsData ()->GetVecArrays ()->Add ("ref2", ref[2]);
    //
    // ********************************************************************* //

    WriteVTKWithCells (store.mesh);
    WriteVTKWithEdges (store.mesh);

    for (Mesh* object : store.listobjects)
    {
        WriteVTKWithCells (object);
        WriteVTKWithEdges (object);
    }

    return EXIT_SUCCESS;
}

