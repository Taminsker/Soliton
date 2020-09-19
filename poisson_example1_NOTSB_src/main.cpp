#include <Soliton>
#include "functions.h"

int main(int argc, char *argv[])
{
    HEADERFUN("main");

    (void)argc;
    (void)argv;

    SolitonRoutines::PrintExecutableName ("POISSON EQUATION - EXECUTABLE 1 - NOT SB");
    SolitonRoutines::PrintMacros ();

    Sto4Sol store;
    InputDatStruct inputdatfile;

    ParseInputDatFile (&inputdatfile, "../input.dat");

    SolitonRoutines::Generation (&store, &inputdatfile);

    // ********************************************************************* //
    BEGIN << COLOR_YELLOW << "Solver" << ENDLINE;

    //    int numPoints = store.mesh->GetNumberOfPoints ();

    //    PlainVector F;
    //    SolitonRoutines::ComputeSecondMember (store.mesh, &F, fun_secondmember);

    //    //    INFOS << F << ENDLINE;

    //    std::vector <Triplet> tripletList;
    //    tripletList.reserve (std::size_t (10 * numPoints));
    //    SolitonRoutines::ComputeGradGrad (store.mesh, &tripletList);

    //    SolitonRoutines::ComputeDirichletNoSBM (store.mesh, &tripletList,  &F, fun_dirichlet, TAG_PHYSICAL::TAG_INLET);
    //    SolitonRoutines::ComputeDirichletNoSBM (store.mesh, &tripletList,  &F, fun_dirichlet, TAG_PHYSICAL::TAG_OUTLET);
    //    SolitonRoutines::ComputeDirichletNoSBM (store.mesh, &tripletList,  &F, fun_dirichlet, TAG_PHYSICAL::TAG_WALL);

    //    //    SolitonRoutines::ComputeDirichletNoSBM (store.mesh, &tripletList,  &F, fun_dirichlet, TAG_PHYSICAL::TAG_DOMAIN, inputdatfile.hsize);

    //    //    SolitonRoutines::ComputeNeumannNoSBM (store.mesh, &tripletList, &F, fun_neumann, TAG_PHYSICAL::TAG_WALL);

    //    INFOS << "Fill the matrix is done !" << ENDLINE;

    //    SparseMatrix A (numPoints, numPoints);
    //    A.setFromTriplets (tripletList.begin (), tripletList.end ());

    //    INFOS <<  "Matrix size : " << A.rows () << ", " << A.cols() << ENDLINE;
    //    std::cout << std::endl;

    //    //    std::ofstream out ("mat.txt");
    //    //    out << A.toDense () << std::endl;
    //    //    out << "\n\n\n\n\n" << std::endl;
    //    //    out << F.transpose ()<< std::endl;
    //    //    out.close ();

    //    PlainVector sol_num;
    //    Solver::AutoDeduceBest (&A, &F, &sol_num);

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

    //           << InitItem({DIRICHLET_BOUNDARY_SBM_ITEM, fun_dirichlet, TAG_INLET, store.listobjects[1]});

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

    PlainVector errorelapercent = GetErrorRelaPercent (store.mesh, &sol_ana, sol_num);
    solstd = PlainVector2Vector (&errorelapercent);
    store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_rela_percent", solstd);

    GetErrorl1 (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);
    GetErrorl2 (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);
    GetErrorlinf (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);
    GetErrorRela (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);

    ENDFUN;

    // ********************************************************************* //

    WriteVTKWithCells (store.mesh, "ex1");

    WriteVTKWithEdges (store.mesh, "ex1");

    for (Mesh* object : store.listobjects)
    {
        WriteVTKWithCells (object, "ex1");
        WriteVTKWithEdges (object, "ex1");
    }

    return EXIT_SUCCESS;
}
