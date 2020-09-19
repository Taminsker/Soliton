#include <Soliton>
#include "functions.h"

int main(int argc, char *argv[])
{
    HEADERFUN("main");

    (void)argc;
    (void)argv;

    SolitonRoutines::PrintExecutableName ("POISSON EQUATION - EXECUTABLE 1 - SB");
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

    UpgradableSolver solver(&store);
    solver.DampingTermOn ();
    solver.TimeSolverOff ();
    
    solver << InitItem({ITEM_SOLVER_TYPE::GRAD_GRAD, _null_fun, TAG_PHYSICAL::TAG_DOMAIN, nullptr, ITEM_SOLVER_DER_TIME::DER_TIME_0})
           << InitItem({ITEM_SOLVER_TYPE::SECOND_MEMBER, fun_secondmember, TAG_PHYSICAL::TAG_DOMAIN, nullptr, ITEM_SOLVER_DER_TIME::DER_TIME_0})
           << InitItem({ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY, fun_dirichlet, TAG_PHYSICAL::TAG_WALL, nullptr, ITEM_SOLVER_DER_TIME::DER_TIME_0})
           << InitItem({ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY, fun_dirichlet, TAG_PHYSICAL::TAG_INLET, nullptr, ITEM_SOLVER_DER_TIME::DER_TIME_0})
           << InitItem({ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY, fun_dirichlet, TAG_PHYSICAL::TAG_OUTLET, nullptr, ITEM_SOLVER_DER_TIME::DER_TIME_0})
           << InitItem({ITEM_SOLVER_TYPE::DIRICHLET_BOUNDARY_SBM, fun_dirichlet, TAG_PHYSICAL::TAG_OBJECT_1, store.listobjects[0], ITEM_SOLVER_DER_TIME::DER_TIME_0});

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

    GetErrorl1 (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);
    GetErrorl2 (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);
    GetErrorlinf (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);
    GetErrorRela (store.mesh, &sol_ana, sol_num, inputdatfile.hsize);

    ENDFUN;
    // ********************************************************************* //

    WriteVTKWithCells (store.mesh, "ex1_sb");

    WriteVTKWithEdges (store.mesh, "ex1_sb");

    for (Mesh* object : store.listobjects)
    {
        WriteVTKWithCells (object, "ex1_sb");
        WriteVTKWithEdges (object, "ex1_sb");
    }

    return EXIT_SUCCESS;
}
