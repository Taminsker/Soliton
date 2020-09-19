#include <Soliton>
#include "functions.h"

int main(int argc, char *argv[])
{
    HEADERFUN("main");

    (void)argc;
    (void)argv;

    SolitonRoutines::PrintExecutableName ("POISSON EQUATION - EXECUTABLE 2 - NOT SB");
    SolitonRoutines::PrintMacros ();

    Sto4Sol store;
    InputDatStruct inputdatfile;

    ParseInputDatFile (&inputdatfile, "../input.dat");

    SolitonRoutines::Generation (&store, &inputdatfile);

    // ********************************************************************* //
    BEGIN << COLOR_YELLOW << "Solver" << ENDLINE;

    int numPoints = store.mesh->GetNumberOfPoints ();

    PlainVector F;
    SolitonRoutines::ComputeSecondMember (store.mesh, &F, f);

    std::vector <Triplet> tripletList;
    tripletList.reserve (std::size_t (10 * numPoints));
    SolitonRoutines::ComputeGradGrad (store.mesh, &tripletList);

    INFOS << "Fill the matrix is done !" << ENDLINE;

    SparseMatrix A (numPoints, numPoints);
    A.setFromTriplets (tripletList.begin (), tripletList.end ());

    INFOS <<  "Matrix size : " << A.rows () << ", " << A.cols() << ENDLINE;
    std::cout << std::endl;


    auto tagvec = store.mesh->GetPointsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec == nullptr)
    {
        INFOS << "no tag physical on points : search under the name " << NAME_TAG_PHYSICAL << "." << BLINKRETURN <<ENDLINE;
        return EXIT_FAILURE;
    }

    std::vector<int> leftright, topbottom;
    leftright.reserve (std::size_t (numPoints));
    topbottom.reserve (std::size_t (numPoints));

    for (int i = 0; i < numPoints; ++i)
    {
        int value = tagvec->vec.at (std::size_t (i));

        if (value == static_cast<int>(TAG_PHYSICAL::TAG_INLET) || value == static_cast<int>(TAG_PHYSICAL::TAG_OUTLET))
            leftright.push_back (i);
        else if (value == static_cast<int>(TAG_PHYSICAL::TAG_WALL))
            topbottom.push_back (i);
    }

    std::sort (topbottom.begin(), topbottom.end());
    topbottom.erase (std::unique (topbottom.begin (), topbottom.end()), topbottom.end());

    std::sort (leftright.begin(), leftright.end());
    leftright.erase (std::unique (leftright.begin (), leftright.end()), leftright.end());

    ImposeDirichlet (store.mesh, &A, &F, g, &topbottom);
    ImposeDirichlet (store.mesh, &A, &F, h, &leftright);

    PlainVector sol_num;
    Solver::AutoDeduceBest (&A, &F, &sol_num);

    std::vector<double> solstd = PlainVector2Vector (&sol_num);

    store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_num", solstd);

    PlainVector sol_ana;
    FunToVec (&sol_ana, store.mesh, u);
    solstd = PlainVector2Vector (&sol_ana);

    store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("sol_ana", solstd);

    PlainVector erroabs = GetErrorAbs (store.mesh, &sol_ana, &sol_num);
    solstd = PlainVector2Vector (&erroabs);
    store.mesh->GetPointsData ()->GetDoubleArrays ()->Add ("err_abs", solstd);

    GetErrorl1 (store.mesh, &sol_ana, &sol_num);
    GetErrorl2 (store.mesh, &sol_ana, &sol_num);
    GetErrorlinf (store.mesh, &sol_ana, &sol_num);
    GetErrorRela (store.mesh, &sol_ana, &sol_num);

    ENDFUN;
    // ********************************************************************* //

    WriteVTKWithCells (store.mesh, "ex2");
    WriteVTKWithEdges (store.mesh, "ex2");

    for (Mesh* object : store.listobjects)
    {
        WriteVTKWithCells (object, "ex2");
        WriteVTKWithEdges (object, "ex2");
    }

    return EXIT_SUCCESS;
}
