#include <Soliton>

double fun_analytic (Point p, double t);
double fun_neumann (Point p, double t);
double fun_dirichlet (Point p, double t);
double fun_secondmember (Point p, double t);

int main(int argc, char *argv[])
{
    HEADERFUN("main");

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
    inputdatfile.fun_neumann = fun_neumann;
    inputdatfile.fun_analytic = fun_analytic;
    inputdatfile.fun_dirichlet = fun_dirichlet;
    inputdatfile.fun_secmember = fun_secondmember;

    ParseInputDatFile (&inputdatfile, argv [1]);

    SolitonRoutines::Generation (&store, &inputdatfile);

    // ********************************************************************* //
    BEGIN << COLOR_YELLOW << "Solver" << ENDLINE;

    auto solver = SolitonRoutines::NewSuperSolver (&store, &inputdatfile);

    solver->SetCoeffPen (inputdatfile.penal);
    solver->SetTime (0.);
    solver->SetTimeStep (inputdatfile.dt);
    solver->SetCollectContribution (inputdatfile.colConcItem);
    solver->Configure ();

    std::vector<std::pair<PlainVector*, INTER>> list_sols_num = solver->Compute (inputdatfile.scheme, true);

    for (std::pair<PlainVector*, INTER> couple : list_sols_num)
    {
        std::vector<double> solstd = PlainVector2Vector (couple.first);

        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_num_"+to_string<INTER>(couple.second), solstd);

        PlainVector sol_ana;
        FunToVec (&sol_ana, store.mesh, fun_analytic);
        solstd = PlainVector2Vector (&sol_ana);

        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_ana_"+to_string<INTER>(couple.second), solstd);

        PlainVector erroabs = GetErrorAbs (store.mesh, &sol_ana, couple.first);
        solstd = PlainVector2Vector (&erroabs);
        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_abs_"+to_string<INTER>(couple.second), solstd);

        PlainVector errorelapercent = GetErrorRelaPercent (store.mesh, &sol_ana, couple.first);
        solstd = PlainVector2Vector (&errorelapercent);
        store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_rela_percent_"+to_string<INTER>(couple.second), solstd);

        GetErrorl1 (store.mesh, &sol_ana, couple.first);
        GetErrorl2 (store.mesh, &sol_ana, couple.first);
        GetErrorlinf (store.mesh, &sol_ana, couple.first);
        GetErrorRela (store.mesh, &sol_ana, couple.first);
    }

    ENDFUN;

    // ********************************************************************* //

    WriteVTKWithCells (store.mesh, "");

    WriteVTKWithEdges (store.mesh, "");

    for (Mesh* object : store.listobjects)
    {
        WriteVTKWithCells (object, "");
        WriteVTKWithEdges (object, "");
    }

    delete solver;

    return EXIT_SUCCESS;
}


double fun_analytic (Point p, double t)
{
    VOID_USE(t);
    return (p.x - p.x * p.x) * (p.y - p.y * p.y);
    //    return 1. * p.x * p.x + 1. * p.y * p.y + 1. * p.x * p.y + 1. * p.x + 1. * p.y + 1.;
    //        return p.x + p.y;
    //    return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
    //    return std::sin(M_PI * p.x) * std::sin(M_PI * p.y);
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

    //        return 0;
    return 2. * (p.x - p.x * p.x + p.y - p.y * p.y);
    //    return -2 * (1. + 1.);
    //    return 2. * M_PI * M_PI * fun_analytic (p, t);
}
