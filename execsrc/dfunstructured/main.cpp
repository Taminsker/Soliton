#include <Soliton>

#include "functions.hpp"

int
main (int argc, char * argv [])
{
    if (argc < 2)
    {
        ERROR << "Use : " << argv [0] << " input.slt" << ENDLINE;
        return EXIT_FAILURE;
    }

    Eigen::setNbThreads (4);
    Eigen::initParallel ();

    SolitonRoutines::PrintExecutableName ("DFUnstructured/Soliton");
    SolitonRoutines::PrintMacros ();

    InputDatStruct inputdatfile;
    ParseInputDatFile (&inputdatfile, argv [1]);

    MeshStorage store;
    SolitonRoutines::Generation (&store, &inputdatfile);

    int                        numPoints    = store.GetMainMesh ()->GetNumberOfPoints ();
    int                        numEdges     = store.GetMainMesh ()->GetNumberOfEdges ();
    PointsData *               ptdata       = store.GetMainMesh ()->GetPointsData ();
    HetContainer<int>::Array * tagInterEdge = store.GetMainMesh ()->GetEdgesData ()->Get<int> (NAME_TAG_INTERSECTION);

    std::vector<real_t> solstd;
    DenseVector     sol_ana, sol_num, err_abs;
    sol_num.setZero (numPoints);
    err_abs.setZero (numPoints);

    FunToVec (&sol_ana, store.GetMainMesh (), &fun_analytic);
    solstd = PlainVector2Vector (&sol_ana);
    ptdata->Add<real_t> ("sol_ana", solstd);

    //** SOLVE PART

    SparseMatrix matrix;

    //* END OF SOLVER

    BEGIN << "Postprocess global" << ENDLINE;
    GetErrorl1 (store.GetMainMesh (), &sol_ana, &sol_num);
    GetErrorl2 (store.GetMainMesh (), &sol_ana, &sol_num);
    GetErrorlinf (store.GetMainMesh (), &sol_ana, &sol_num);
    GetErrorRela (store.GetMainMesh (), &sol_ana, &sol_num);

    solstd = PlainVector2Vector (&sol_num);
    ptdata->Add<real_t> ("sol_num", solstd);
    solstd = PlainVector2Vector (&err_abs);
    ptdata->Add<real_t> ("err_abs", solstd);

    ENDFUN;

    // Write part
    WriteVTKWithCells (store.GetMainMesh ());
    WriteVTKWithEdges (store.GetMainMesh ());

    for (Mesh * object : *store.GetListOfObjects ())
    {
        WriteVTKWithCells (object);
        WriteVTKWithEdges (object);
    }

    return EXIT_SUCCESS;
}
