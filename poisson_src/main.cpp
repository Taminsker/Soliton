
#include <iostream>
#include <variant>
#include <Soliton>

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    auto _err_ = system ("clear");
    (void)_err_;

    std::cout << REVERSE << COLOR_BLUE << "--------------------------------------------------------" << ENDLINE;
    std::cout << REVERSE << COLOR_BLUE << "                     POISSON EQUATION                   " << ENDLINE;
    std::cout << REVERSE << COLOR_BLUE << "--------------------------------------------------------" << ENDLINE << ENDLINE;
#ifdef VERBOSE
    INFOS << "Verbose is ON" << ENDLINE;
#else
    INFOS << "Verbose is OFF" << ENDLINE;
#endif

#ifdef DEBUG
    INFOS << "Debug is ON" << ENDLINE;
#else
    INFOS << "DEBUG is OFF" << ENDLINE;
#endif
    std::cout << ENDLINE;

    HEADERFUN("main");

    Sto4Sol store;
    InputDatStruct inputdatfile;
    SOLITON_RETURN error;

    error = ParseInputDatFile (&inputdatfile, "../input.dat");
    USE_SOLITON_RETURN(error);

    error = ParseMSH (store.mesh, GMSH::Generate (&inputdatfile), false);
    USE_SOLITON_RETURN(error);

    error = ObjectGenerator (&inputdatfile, &store);
    USE_SOLITON_RETURN(error);

    //    size_t numObjects = store.listobjects.size ();
    //    store.mesh->Print();

    for (Mesh* object : store.listobjects)
    {
        //        object->Print ();

        error = AddLevelSetAndTag (store.mesh, object);
        USE_SOLITON_RETURN(error);

        error = AddLevelSetAndTag (object, store.mesh);
        USE_SOLITON_RETURN(error);
    }

    error = BuildSurrogateDomains (&store);
    USE_SOLITON_RETURN(error);

    //    store.mesh->Print();

    // ********************************************************************* //
    BEGIN << COLOR_YELLOW << "Solver Poisson" << ENDLINE;

    DFStruct dfstruct;
    FEStruct festruct;

    PlainVector F;
    FunToVec (&F, store.mesh, 0.0);

    int numPoints = store.mesh->GetNumberOfPoints ();
    int numCells = store.mesh->GetNumberOfCells ();

    std::vector <Eigen::Triplet<double> > tripletList;
    tripletList.reserve (std::size_t (2 * numPoints));

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = store.mesh->GetCell (cellId);
//        constexpr VTKCellType type = VTK_LINE;
        const VTKCellType type = cell->GetTypeVTK ();


        auto* element = festruct.GetElementFor (type);

//        tripletList.push_back(T(i,j,v_ij));
    }

    SparseMatrix A (numPoints, numPoints);
    A.setFromTriplets (tripletList.begin (), tripletList.end ());



    ENDFUN;
    // ********************************************************************* //

    error = WriterVTK::WithCells (store.mesh);
    USE_SOLITON_RETURN(error);

    error = WriterVTK::WithEdges (store.mesh);
    USE_SOLITON_RETURN(error);


    for (Mesh* object : store.listobjects)
    {
        error = WriterVTK::WithCells (object);
        USE_SOLITON_RETURN(error);
        error = WriterVTK::WithEdges (object);
        USE_SOLITON_RETURN(error);
    }




    return SOLITON_SUCCESS;
}
