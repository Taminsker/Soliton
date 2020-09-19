#include "routines.h"

void SolitonRoutines::PrintMacros ()
{
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

//#ifdef PRINTHASHMAP
//    INFOS << "PRINTHASHMAP is ON" << ENDLINE;
//#else
//    INFOS << "PRINTHASHMAP is OFF" << ENDLINE;
//#endif

    std::cout << ENDLINE;

    return;
}

void SolitonRoutines::PrintExecutableName (std::string executablename)
{
    std::size_t n = 20;
    std::size_t s = executablename.size ();

    std::cout << REVERSE << COLOR_BLUE;

    // first line
    for (std::size_t i = 0; i < 2*n+s; ++i)
        std::cout << "-";
    std::cout << "\n" << std::flush;

    // Name
    for (std::size_t i = 0; i < n; ++i)
        std::cout << " ";

    std::cout << executablename;
    for (std::size_t i = 0; i < n; ++i)
        std::cout << " ";

    std::cout << "\n" << std::flush;


    // second line
    for (std::size_t i = 0; i < 2*n+s; ++i)
        std::cout << "-";
    std::cout << "\n" << ENDLINE;

    return;
}

void SolitonRoutines::Generation (Sto4Sol* store, InputDatStruct* inputdatfile)
{

    if (inputdatfile->gen)
    {
        ParseMSH (store->mesh, GenerateWithGMSH (inputdatfile), false);
        store->mesh->SetPrescribedSize (inputdatfile->hsize);
    }
    else
        ParseMSH (store->mesh, inputdatfile->filename_msh, true);


    /*************************/

    TransfertTagPhysical (store->mesh);

    if (inputdatfile->damping)
        ComputeDampingArea (store->mesh, TAG_PHYSICAL::TAG_OUTLET, inputdatfile->hsize);


    ObjectGenerator (inputdatfile, store);

#ifdef VERBOSE
    STATUS << "number of objects : " << store->listobjects.size () << ENDLINE;
    ENDFUN;
#endif

    for (std::size_t objectId = 0; objectId < store->listobjects.size (); ++objectId)
    {
        Mesh* object = store->listobjects.at (objectId);

#ifdef VERBOSE
        object->Print ();
#endif
        // MESH -> OBJECT
        AddLevelSetBetween (store->mesh, object);
        TagCellsFromLevelSet (store->mesh, object);
        TagEdgesFromTagCells (store->mesh, object);

#ifndef VERBOSE
        ENDFUN;
#endif

        // OBJECT -> MESH
        AddLevelSetBetween (object, store->mesh);
        TagCellsFromLevelSet (object, store->mesh);
        TagEdgesFromTagCells (object, store->mesh);

#ifndef VERBOSE
        ENDFUN;
#endif
        // Displacement vector
        BuildDisplacementVectorsBounds (store->mesh, object);
    }
#ifdef VERBOSE
    store->mesh->Print();
#endif

    return;
}

void SolitonRoutines::ComputeGradGrad (Mesh* mesh, std::vector<Triplet>* listTriplets)
{
    // * <grad(phi_i), grad(phi_j)>
    FEStore festore;
    QuadStore quadstore;

    int numCells = mesh->GetNumberOfCells ();

    FELocalInfos loc;

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);

        if (cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_LINE || cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_QUADRATIC_EDGE)
            continue;

        std::vector <Point *>* listpts = cell->GetPoints ();
        std::size_t numCellPoints = listpts->size ();
        FEBase* element = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);
        QuadStore::QuadObject* quad = quadstore.Get (element);

        for (std::size_t i = 0; i < numCellPoints; ++i)
        {
            Point* pt_i = listpts->at (i);
            int pt_iIdx = pt_i->GetGlobalIndex ();

            for (std::size_t j = 0; j < numCellPoints; ++j)
            {
                Point* pt_j = listpts->at (j);
                int pt_jIdx = pt_j->GetGlobalIndex ();

                double coeff = 0.;
                Point v1, v2;

                std::function <Point(Point*)> grad_phi_i = element->GetGradPhi (i);
                std::function <Point(Point*)> grad_phi_j = element->GetGradPhi (j);

                for (std::size_t k = 0; k < quad->npts; ++k)
                {
                    element->LocalCompute (&quad->pts [k], &loc);

                    v1 = loc.JacInvT * grad_phi_i (&quad->pts [k]);
                    v2 = loc.JacInvT * grad_phi_j (&quad->pts [k]);

                    coeff += loc.detJac * quad->w [k] * (v1 | v2);
                }

                listTriplets->push_back({pt_iIdx, pt_jIdx, coeff});
            }
        }
    }

    return;
}

void SolitonRoutines::ComputeSecondMember (Mesh* mesh, PlainVector* vec, double (*f)(Point, double))
{
    // * <f, phi_j>

    vec->setZero (mesh->GetNumberOfPoints ());

    FEStore festore;
    QuadStore quadstore;

    int numCells = mesh->GetNumberOfCells ();
    FELocalInfos loc;

    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);

        if (cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_LINE || cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_QUADRATIC_EDGE)
            continue;

        std::vector <Point *>* listpts = cell->GetPoints ();
        std::size_t numCellPoints = listpts->size ();
        FEBase* element = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);
        QuadStore::QuadObject* quad = quadstore.Get (element);

        for (std::size_t j = 0; j < numCellPoints; ++j)
        {
            Point* pt_j = listpts->at (j);
            int pt_jIdx = pt_j->GetGlobalIndex ();
            double coeff = 0;
            std::function <double(Point*)> phi_j = element->GetPhi (j);

            for (std::size_t k = 0; k < quad->npts; ++k)
            {
                element->LocalCompute (&quad->pts [k], &loc);
                Point npt = element->TransformRefToEle (&quad->pts [k]);

                coeff += loc.detJac * quad->w [k] * f(npt, 0) * phi_j (&quad->pts [k]);
            }

            vec->coeffRef (pt_jIdx) += coeff;
        }
    }

    return;
}



void SolitonRoutines::ComputeNeumannNoSBM (Mesh* mesh, std::vector<Triplet>* listTriplets, PlainVector* F, double (*f)(Point, double), TAG_PHYSICAL tagToApply, double hPerp, double t)
{
    // * +<phi_j, grad(u).n - t_n>

    // decouple
    // * A : +<phi_j, grad(u).n>
    // * B : +<phi_j, t_n>

    (void)hPerp;

    FEStore festore;
    QuadStore quadstore;

    INFOS << "Impose Neumann : tag " << tagToApply << ENDLINE;
    int numCells = mesh->GetNumberOfCells ();

    auto tagedgevec = mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagedgevec == nullptr)
    {
        ERROR << "no array with the name : " << NAME_TAG_PHYSICAL << "." << BLINKRETURN << ENDLINE;
        return;
    }

    auto edgeNormals = mesh->GetEdgesData ()->GetVecArrays ()->Get (NAME_NORMAL_ON_EDGES);

    if (edgeNormals == nullptr)
    {
        ERROR << "no array with the name : " << NAME_NORMAL_ON_EDGES << "." << BLINKRETURN << ENDLINE;
        return;
    }

    FELocalInfos locEdge, locCell;


    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);

        if (cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_LINE)
            continue;

        std::vector<Edge*>* listEdgesOnCell = cell->GetEdges ();
        std::size_t numOfEdgesOnCell = listEdgesOnCell->size ();
        std::vector<Point*>* listPointsOnCell = cell->GetPoints ();
        std::size_t numOfPointsOnCell = listPointsOnCell->size ();

        FEBase* FEElement4Cell = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);

        for (std::size_t i = 0; i < numOfPointsOnCell; ++i)
        {
            Point* pt_i = listPointsOnCell->at (i);
            int idGlobal_i = pt_i->GetGlobalIndex ();
            LambdaOnPoint<Point> grad_phi_i = FEElement4Cell->GetGradPhi (i);
            LambdaOnPoint<double> phi_i = FEElement4Cell->GetPhi (i);


            for (std::size_t j = 0; j < numOfPointsOnCell; ++j)
            {
                Point* pt_j = listPointsOnCell->at (j);
                int idGlobal_j = pt_j->GetGlobalIndex ();
                LambdaOnPoint<Point> grad_phi_j = FEElement4Cell->GetGradPhi (j);
                LambdaOnPoint<double> phi_j = FEElement4Cell->GetPhi (j);


                for (std::size_t edgeLocId = 0; edgeLocId < numOfEdgesOnCell; ++edgeLocId)
                {
                    Edge* edge = listEdgesOnCell->at (edgeLocId);
                    int idGlobal_edge = edge->GetGlobalIndex ();


                    if (tagedgevec->vec.at (std::size_t (idGlobal_edge)) != static_cast<int>(tagToApply))
                        continue;

                    Edge* edgeOnCellRef = FEElement4Cell->GetEdge (edgeLocId);

                    Point* normal = edgeNormals->vec.at (std::size_t (idGlobal_edge));
                    double prodscal = ((*edge->GetCentroid () - *cell->GetCentroid ()) | *normal);

                    if (prodscal < 0)
                        *normal = -1 * *normal;

                    FEBase* FEElement4Edge = festore.GetElementFor (edgeOnCellRef);
                    QuadStore::QuadObject* quad4Edge = quadstore.Get (FEElement4Edge);

                    double coeff1, coeff2;
                    Point v1, v2;
                    double value1, value2;
                    double pdetJac;
                    Matrix3x3 pJacInvT;

                    coeff1 = 0.;
                    coeff2 = 0.;

                    for (std::size_t k = 0; k < quad4Edge->npts; ++k)
                    {
                        double w_k = quad4Edge->w [k];
                        //                        Point pt_int_k = quad4Edge->pts [k];

                        FEElement4Edge->LocalCompute (&quad4Edge->pts [k], &locEdge);

                        Point pt_k = FEElement4Edge->TransformRefToEle (&quad4Edge->pts [k]);

                        FEElement4Cell->LocalCompute (&pt_k, &locCell);

                        pdetJac = locEdge.detJac * locCell.detJac;
                        pJacInvT = locEdge.JacInvT * locCell.JacInvT;

                        // * A : +<phi_j, grad(phi_i).n>

                        value1 = phi_j(&pt_k);
                        value2 = ((pJacInvT * grad_phi_i(&pt_k)) | (pJacInvT * *normal));
                        coeff1 += pdetJac * w_k * (value1 * value2);


                        if (i == 0)
                        {
                            // Second membre

                            // * : +<phi_j, t_n>

                            value1 = phi_j (&pt_k);
                            value2 = f(FEElement4Cell->TransformRefToEle (&pt_k), t);
                            coeff2 += pdetJac * w_k * (value1 * value2);
                        }
                    }

                    listTriplets->push_back({idGlobal_i, idGlobal_j, coeff1});

                    if (i == 0)
                        F->coeffRef (idGlobal_j) += coeff2;
                }
            }
        }
    }
    return;
}


void SolitonRoutines::ComputeDirichletNoSBM (Mesh* mesh, std::vector<Triplet>* listTriplets, PlainVector* F, double (*f)(Point, double), TAG_PHYSICAL tagToApply, double hsz, double t)
{

    // * -<phi_j, grad(u).n>
    // * - <grad(phi_j).n, u - u_D>
    // * + <alpha * phi_j, u - u_D>

    // decouple
    // * A : -<phi_j, grad(u).n>
    // * A : - <grad(phi_j).n, u>
    // * A : + <alpha * phi_j, u>
    //
    // * B : - <grad(phi_j).n, u_D>
    // * B : + <alpha * phi_j, u_D>

    FEStore festore;
    QuadStore quadstore;

    double max_aij = listTriplets->front ().value ();
    for (Triplet triplet : *listTriplets)
        max_aij = std::max (max_aij, triplet.value ());

    double alpha = 20 / std::abs(hsz);

    INFOS << "Impose Dirichlet : tag " << tagToApply << ", penal. coeff. = " << alpha << ENDLINE;
    int numCells = mesh->GetNumberOfCells ();

    auto tagedgevec = mesh->GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagedgevec == nullptr)
    {
        ERROR << "no array with the name : " << NAME_TAG_PHYSICAL << "." << BLINKRETURN << ENDLINE;
        return;
    }

    FELocalInfos locEdge, locCell;


    for (int cellId = 0; cellId < numCells; ++cellId)
    {
        Cell* cell = mesh->GetCell (cellId);

        if (cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_LINE || cell->GetTypeVTK () == VTK_CELL_TYPE::VTK_VERTEX)
            continue;

        std::vector<Edge*>* listEdgesOnCell = cell->GetEdges ();
        std::size_t numOfEdgesOnCell = listEdgesOnCell->size ();
        std::vector<Point*>* listPointsOnCell = cell->GetPoints ();
        std::size_t numOfPointsOnCell = listPointsOnCell->size ();

        FEBase* FEElement4Cell = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);

        for (std::size_t i = 0; i < numOfPointsOnCell; ++i)
        {
            Point* pt_i = listPointsOnCell->at (i);
            int idGlobal_i = pt_i->GetGlobalIndex ();
            LambdaOnPoint<Point> grad_phi_i = FEElement4Cell->GetGradPhi (i);
            LambdaOnPoint<double> phi_i = FEElement4Cell->GetPhi (i);

            for (std::size_t j = 0; j < numOfPointsOnCell; ++j)
            {
                Point* pt_j = listPointsOnCell->at (j);
                int idGlobal_j = pt_j->GetGlobalIndex ();
                LambdaOnPoint<Point> grad_phi_j = FEElement4Cell->GetGradPhi (j);
                LambdaOnPoint<double> phi_j = FEElement4Cell->GetPhi (j);


                for (std::size_t edgeLocId = 0; edgeLocId < numOfEdgesOnCell; ++edgeLocId)
                {
                    Edge* edge = listEdgesOnCell->at (edgeLocId);
                    int idGlobal_edge = edge->GetGlobalIndex ();


                    if (tagedgevec->vec.at (std::size_t (idGlobal_edge)) != static_cast<int>(tagToApply))
                        continue;

                    Edge* edgeOnCellRef = FEElement4Cell->GetEdge (edgeLocId);


                    FEBase* FEElement4Edge = festore.GetElementFor (edgeOnCellRef);
                    QuadStore::QuadObject* quad4Edge = quadstore.Get (FEElement4Edge);


                    double coeff1, coeff2, ck;
                    Point v1, v2;
                    Matrix3x3 pJacInvT;
                    double phi_i_eval, phi_j_eval;
                    Point grad_phi_i_eval, grad_phi_j_eval, normal_eval;
                    double u_D;

                    coeff1 = 0.;
                    coeff2 = 0.;
                    ck = 0.;
                    phi_i_eval = 0.;
                    phi_j_eval = 0.;
                    grad_phi_i_eval = {0., 0., 0.};
                    grad_phi_j_eval = {0., 0., 0.};
                    normal_eval = {0., 0., 0.};
                    u_D = 0.;

                    for (std::size_t k = 0; k < quad4Edge->npts; ++k)
                    {

                        Point pt_k_oncell = FEElement4Edge->TransformRefToEle (&quad4Edge->pts [k]);
                        Point pt_k_real = FEElement4Cell->TransformRefToEle (&pt_k_oncell);

                        FEElement4Edge->LocalCompute (&quad4Edge->pts [k], &locEdge);
                        FEElement4Cell->LocalCompute (&pt_k_oncell, &locCell);

                        ck = locEdge.detJac * locCell.detJac * quad4Edge->w [k];

                        phi_i_eval = phi_i (&pt_k_oncell);
                        phi_j_eval = phi_j (&pt_k_oncell);
                        grad_phi_i_eval = locEdge.JacInvT * locCell.JacInvT * grad_phi_i (&pt_k_oncell);
                        grad_phi_j_eval = locEdge.JacInvT * locCell.JacInvT * grad_phi_j (&pt_k_oncell);
                        //                        normal_eval = JacInvTEdge * JacInvTCell * *normal;
                        normal_eval = {0, 1, 0};

                        u_D = f(pt_k_real, t);

                        // * A : -<phi_j, grad(phi_i).n>

                        coeff1 -= ck * phi_j_eval * (grad_phi_i_eval | normal_eval);

                        // * A : - <grad(phi_j).n, phi_i>

                        coeff1 -= ck * (grad_phi_j_eval | normal_eval) * phi_i_eval;

                        // * A : + <alpha * phi_j, phi_i>

                        coeff1 += ck * alpha * phi_j_eval * phi_i_eval;

                        if (i == 0)
                        {
                            // Second membre

                            // * - <grad(phi_j).n, u_D>
                            coeff2 -= ck * (grad_phi_j_eval | normal_eval) * u_D;

                            // * + <alpha * phi_j, u_D>
                            coeff2 += ck * alpha * phi_j_eval * u_D;
                        }
                    }

                    listTriplets->push_back({idGlobal_i, idGlobal_j, coeff1});

                    if (i == 0)
                        F->coeffRef (idGlobal_j) += coeff2;
                }
            }
        }
    }
    return;
}



//void SolitonRoutines::ComputeDirichletSBM (Mesh* mesh, Mesh* object, std::vector<Triplet>* listTriplets, PlainVector* F, double (*f)(Point, double), double hsz, double t)
//{

//    // * -<phi_j, grad(u).n>
//    // * - <grad(phi_j).n, u - u_D>
//    // * + <alpha * phi_j, u - u_D>

//    // decouple
//    // * A : -<phi_j, grad(u).n>
//    // * A : - <grad(phi_j).n, u>
//    // * A : + <alpha * phi_j, u>
//    //
//    // * B : - <grad(phi_j).n, u_D>
//    // * B : + <alpha * phi_j, u_D>

//    FEStore festore;
//    QuadStore quadstore;

//    double max_aij = listTriplets->front ().value ();
//    for (Triplet triplet : *listTriplets)
//        max_aij = std::max (max_aij, triplet.value ());

//    double alpha = 20 / std::abs(hsz);

//    INFOS << "Impose Dirichlet : name of vector of tag surrogate " << object->GetName ()+  NAME_TAG_SURROGATE << ", penal. coeff. = " << alpha << ENDLINE;
//    int numCells = mesh->GetNumberOfCells ();

//    HetInt::Array* tagedgevec = mesh->GetEdgesData ()->GetIntArrays ()->Get (object->GetName ()+  NAME_TAG_SURROGATE);

//    if (tagedgevec == nullptr)
//    {
//        ERROR << "no array with the name : " << object->GetName ()+  NAME_TAG_SURROGATE << "." << BLINKRETURN << ENDLINE;
//        return;
//    }

//    HetPointptr::Array* dvec = mesh->GetPointsData ()->GetVecArrays ()->Get (object->GetName () + NAME_DISPLACEMENT_VECTOR);
//    if (dvec == nullptr)
//    {
//        ERROR << "no array with the name : " << object->GetName () + NAME_DISPLACEMENT_VECTOR << "." << BLINKRETURN << ENDLINE;
//        return;
//    }

//    FELocalInfos locEdge, locCell;


//    for (int cellId = 0; cellId < numCells; ++cellId)
//    {
//        Cell* cell = mesh->GetCell (cellId);

//        if (cell->GetTypeVTK () == VTK_LINE || cell->GetTypeVTK () == VTK_VERTEX)
//            continue;

//        std::vector<Edge*>* listEdgesOnCell = cell->GetEdges ();
//        std::size_t numOfEdgesOnCell = listEdgesOnCell->size ();
//        std::vector<Point*>* listPointsOnCell = cell->GetPoints ();
//        std::size_t numOfPointsOnCell = listPointsOnCell->size ();

//        FEBase* FEElement4Cell = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);

//        for (std::size_t i = 0; i < numOfPointsOnCell; ++i)
//        {
//            Point* pt_i = listPointsOnCell->at (i);
//            int idGlobal_i = pt_i->GetGlobalIndex ();
//            LambdaOnPoint<Point> grad_phi_i = FEElement4Cell->gradPhi [i];
//            LambdaOnPoint<double> phi_i = FEElement4Cell->Phi [i];

//            for (std::size_t j = 0; j < numOfPointsOnCell; ++j)
//            {
//                Point* pt_j = listPointsOnCell->at (j);
//                int idGlobal_j = pt_j->GetGlobalIndex ();
//                LambdaOnPoint<Point> grad_phi_j = FEElement4Cell->gradPhi [j];
//                LambdaOnPoint<double> phi_j = FEElement4Cell->Phi [j];


//                for (std::size_t edgeLocId = 0; edgeLocId < numOfEdgesOnCell; ++edgeLocId)
//                {
//                    Edge* edge = listEdgesOnCell->at (edgeLocId);
//                    int idGlobal_edge = edge->GetGlobalIndex ();


//                    if (tagedgevec->vec.at (std::size_t (idGlobal_edge)) != tagToApply)
//                        continue;

//                    Edge* edgeOnCellRef = FEElement4Cell->edges [edgeLocId];

////                    Point* normal = edgeNormals->vec.at (std::size_t (idGlobal_edge));
////                    double prodscal = ((*edge->GetCentroid () - *cell->GetCentroid ()) | *normal);

////                    if (prodscal < 0)
////                        *normal = -1 * *normal;

//                    FEBase* FEElement4Edge = festore.GetElementFor (edgeOnCellRef);
//                    QuadStore::QuadObject* quad4Edge = quadstore.Get (FEElement4Edge);


//                    double coeff1, coeff2, ck;
//                    Point v1, v2;
//                    Matrix3x3 pJacInvT;
//                    double phi_i_eval, phi_j_eval;
//                    Point grad_phi_i_eval, grad_phi_j_eval, normal_eval;
//                    double u_D;

//                    coeff1 = 0.;
//                    coeff2 = 0.;
//                    ck = 0.;
//                    phi_i_eval = 0.;
//                    phi_j_eval = 0.;
//                    grad_phi_i_eval = {0., 0., 0.};
//                    grad_phi_j_eval = {0., 0., 0.};
//                    normal_eval = {0., 0., 0.};
//                    u_D = 0.;

//                    for (std::size_t k = 0; k < quad4Edge->npts; ++k)
//                    {

//                        Point pt_k_oncell = FEElement4Edge->TransformRefToEle (&quad4Edge->pts [k]);
//                        Point pt_k_real = FEElement4Cell->TransformRefToEle (&pt_k_oncell);

//                        FEElement4Edge->LocalCompute (&quad4Edge->pts [k], &locEdge);
//                        FEElement4Cell->LocalCompute (&pt_k_oncell, &locCell);

//                        ck = locEdge.detJac * locCell.detJac * quad4Edge->w [k];

//                        phi_i_eval = phi_i (&pt_k_oncell);
//                        phi_j_eval = phi_j (&pt_k_oncell);
//                        grad_phi_i_eval = locEdge.JacInvT * locCell.JacInvT * grad_phi_i (&pt_k_oncell);
//                        grad_phi_j_eval = locEdge.JacInvT * locCell.JacInvT * grad_phi_j (&pt_k_oncell);
//                        //                        normal_eval = JacInvTEdge * JacInvTCell * *normal;
//                        normal_eval = {0, 1, 0};

//                        u_D = f(pt_k_real, t);

//                        // * A : -<phi_j, grad(phi_i).n>

//                        coeff1 -= ck * phi_j_eval * (grad_phi_i_eval | normal_eval);

//                        // * A : - <grad(phi_j).n, phi_i>

//                        coeff1 -= ck * (grad_phi_j_eval | normal_eval) * phi_i_eval;

//                        // * A : + <alpha * phi_j, phi_i>

//                        coeff1 += ck * alpha * phi_j_eval * phi_i_eval;

//                        if (i == 0)
//                        {
//                            // Second membre

//                            // * - <grad(phi_j).n, u_D>
//                            coeff2 -= ck * (grad_phi_j_eval | normal_eval) * u_D;

//                            // * + <alpha * phi_j, u_D>
//                            coeff2 += ck * alpha * phi_j_eval * u_D;
//                        }

////                        if (i == 0 && idGlobal_i == 283)
////                            INFOS << "(i) 283 ck = " << coeff1  << " " << cellId << ENDLINE;
////                        if (j == 0 && idGlobal_j == 283)
////                            INFOS << "(j) 283 ck = " << coeff1  << " " << cellId << ENDLINE;

////                        if (i == 0 && idGlobal_i == 279)
////                            INFOS << "(i) 279 ck = " << coeff1 << " " << cellId << ENDLINE;
////                        if (j == 0 && idGlobal_j == 279)
////                            INFOS << "(j) 279 ck = " << coeff1 << " " << cellId << ENDLINE;
//                    }

//                    listTriplets->push_back({idGlobal_i, idGlobal_j, coeff1});

//                    if (i == 0)
//                        F->coeffRef (idGlobal_j) += coeff2;
//                }
//            }
//        }
//    }
//    return;
//}
