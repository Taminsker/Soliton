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
    // INFOS << "PRINTHASHMAP is ON" << ENDLINE;
    //#else
    // INFOS << "PRINTHASHMAP is OFF" << ENDLINE;
    //#endif

    std::cout << ENDLINE;

    return;
}

void SolitonRoutines::PrintExecutableName (std::string executablename)
{
    ul_t n = 20;
    ul_t s = executablename.size ();

    std::cout << REVERSE << COLOR_BLUE;

    // first line
    for (ul_t i = 0; i < 2*n+s; ++i)
        std::cout << "-";
    std::cout << "\n" << std::flush;

    // Name
    for (ul_t i = 0; i < n; ++i)
        std::cout << " ";

    std::cout << executablename;
    for (ul_t i = 0; i < n; ++i)
        std::cout << " ";

    std::cout << "\n" << std::flush;


    // second line
    for (ul_t i = 0; i < 2*n+s; ++i)
        std::cout << "-";
    std::cout << "\n" << ENDLINE;

    return;
}

void SolitonRoutines::Generation (MeshStorage* store, InputDatStruct* inputdatfile)
{

    if (inputdatfile->filename_msh.empty ())
    {
        ParseMSH (store->GetMainMesh (), GenerateWithGMSH (inputdatfile), false);
        store->GetMainMesh ()->SetPrescribedSize (inputdatfile->hsize);
    }
    else
        ParseMSH (store->GetMainMesh (), inputdatfile->filename_msh, true);

    /*************************/

    ComputeTagPhysical (store->GetMainMesh (), inputdatfile);

    if (inputdatfile->damping)
        ComputeDampingArea (store->GetMainMesh (), PHYS::OUTLET, inputdatfile->hsize);


    ObjectGenerator (inputdatfile, store);

#ifdef VERBOSE
    STATUS << "number of objects : " << store->GetListOfObjects ()->size () << ENDLINE;
    ENDFUN;
#endif

    // Definitions
    ul_t numCells = static_cast<ul_t>(store->GetMainMesh ()->GetNumberOfCells ());
    ul_t numEdges = static_cast<ul_t>(store->GetMainMesh ()->GetNumberOfEdges ());
    ul_t numPoints = static_cast<ul_t>(store->GetMainMesh ()->GetNumberOfPoints ());

    // Tags definitions
    std::vector<int> interTagsCells, interTagsEdges, interTagsPoints;
    interTagsCells.resize (numCells, int(INTER::DEFAULT));
    interTagsEdges.resize (numEdges, int(INTER::DEFAULT));
    interTagsPoints.resize (numPoints, int(INTER::DEFAULT));

    std::vector<Point*> displacementVector (numPoints);
    for (ul_t i = 0; i < numPoints; ++i)
        displacementVector.at (i) = new Point ();

    for (ul_t objectId = 0; objectId < store->GetListOfObjects ()->size (); ++objectId)
    {
        Mesh* object = store->GetMeshObjectAt (objectId);

#ifdef VERBOSE
        object->Print ();
#endif
        // MESH -> OBJECT
        AddLevelSetBetween (store->GetMainMesh (), object);
        TagCellsFromLevelSet (store->GetMainMesh (), object);
        TagEdgesFromTagCells (store->GetMainMesh (), object);

#ifndef VERBOSE
        ENDFUN;
#endif

        // OBJECT -> MESH
        AddLevelSetBetween (object, store->GetMainMesh ());
        TagCellsFromLevelSet (object, store->GetMainMesh ());
        TagEdgesFromTagCells (object, store->GetMainMesh ());

#ifndef VERBOSE
        ENDFUN;
#endif
        // Displacement vector
        BuildDisplacementVectorsBounds (store->GetMainMesh (), object);

        // Reduce tags
        HetContainer<int>::Array* currentTagsInterCells = store->GetMainMesh ()->GetCellsData ()->Get<int> (object->GetName () + NAME_TAG_INTERSECTION);
        HetContainer<int>::Array* currentTagsInterEdges = store->GetMainMesh ()->GetEdgesData ()->Get<int> (object->GetName () + NAME_TAG_INTERSECTION);
        HetContainer<Point*>::Array* currentDisplacementPoints = store->GetMainMesh ()->GetPointsData ()->Get<Point*> (object->GetName () + NAME_DISPLACEMENT_VECTOR);

        for (ul_t i = 0; i < numCells; ++i)
            interTagsCells.at (i) |= currentTagsInterCells->vec.at (i);

        for (ul_t i = 0; i < numEdges; ++i)
            interTagsEdges.at (i) |= currentTagsInterEdges->vec.at (i);

        for (ul_t i = 0; i < numPoints; ++i)
            *displacementVector.at (i) += *currentDisplacementPoints->vec.at (i);

    }

    for (int &value : interTagsCells)
        if (value == int(INTER::MIXED))
            value = int(INTER::IN);

    store->GetMainMesh ()->GetCellsData ()->Add<int> (NAME_TAG_INTERSECTION, interTagsCells);
    store->GetMainMesh ()->GetEdgesData ()->Add<int> (NAME_TAG_INTERSECTION, interTagsEdges);
    // store->GetMainMesh ()->GetPointsData ()->Add<int> (NAME_TAG_INTERSECTION, interTagsPoints);
    store->GetMainMesh ()->GetPointsData ()->Add<Point*> (NAME_DISPLACEMENT_VECTOR, displacementVector);

#ifdef VERBOSE
    store->GetMainMesh ()->Print();
#endif

    return;
}

void SolitonRoutines::ComputeGradGrad (Mesh* mesh, std::vector<Triplet_eig>* listTriplets)
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
        ul_t numCellPoints = listpts->size ();
        FEBase* element = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);
        QuadStore::QuadObject* quad = quadstore.Get (element);

        for (ul_t i = 0; i < numCellPoints; ++i)
        {
            Point* pt_i = listpts->at (i);
            int pt_iIdx = pt_i->GetGlobalIndex ();

            for (ul_t j = 0; j < numCellPoints; ++j)
            {
                Point* pt_j = listpts->at (j);
                int pt_jIdx = pt_j->GetGlobalIndex ();

                real_t coeff = 0.;
                Point v1, v2;

                std::function <Point(Point*)> grad_phi_i = element->GetGradPhi (i);
                std::function <Point(Point*)> grad_phi_j = element->GetGradPhi (j);

                for (ul_t k = 0; k < quad->npts; ++k)
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

void SolitonRoutines::ComputeSecondMember (Mesh* mesh, PlainVector_eig *vec, real_t (*f)(Point, real_t))
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
        ul_t numCellPoints = listpts->size ();
        FEBase* element = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);
        QuadStore::QuadObject* quad = quadstore.Get (element);

        for (ul_t j = 0; j < numCellPoints; ++j)
        {
            Point* pt_j = listpts->at (j);
            int pt_jIdx = pt_j->GetGlobalIndex ();
            real_t coeff = 0;
            std::function <real_t(Point*)> phi_j = element->GetPhi (j);

            for (ul_t k = 0; k < quad->npts; ++k)
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



void SolitonRoutines::ComputeNeumannNoSBM (Mesh* mesh, std::vector<Triplet_eig>* listTriplets, PlainVector_eig* F, real_t (*f)(Point, real_t), PHYS tagToApply, real_t, real_t t)
{
    // * +<phi_j, grad(u).n - t_n>

    // decouple
    // * A : +<phi_j, grad(u).n>
    // * B : +<phi_j, t_n>

    FEStore festore;
    QuadStore quadstore;

    INFOS << "Impose Neumann : tag " << to_string(tagToApply) << ENDLINE;
    int numCells = mesh->GetNumberOfCells ();

    auto tagedgevec = mesh->GetEdgesData ()->Get<int> (NAME_TAG_PHYSICAL);

    if (tagedgevec == nullptr)
    {
        ERROR << "no array with the name : " << NAME_TAG_PHYSICAL << "." << BLINKRETURN << ENDLINE;
        return;
    }

    auto edgeNormals = mesh->GetEdgesData ()->Get<Point*> (NAME_NORMAL_ON_EDGES);

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
        ul_t numOfEdgesOnCell = listEdgesOnCell->size ();
        std::vector<Point*>* listPointsOnCell = cell->GetPoints ();
        ul_t numOfPointsOnCell = listPointsOnCell->size ();

        FEBase* FEElement4Cell = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);

        for (ul_t i = 0; i < numOfPointsOnCell; ++i)
        {
            Point* pt_i = listPointsOnCell->at (i);
            int idGlobal_i = pt_i->GetGlobalIndex ();
            LambdaOnPoint<Point> grad_phi_i = FEElement4Cell->GetGradPhi (i);
            LambdaOnPoint<real_t> phi_i = FEElement4Cell->GetPhi (i);


            for (ul_t j = 0; j < numOfPointsOnCell; ++j)
            {
                Point* pt_j = listPointsOnCell->at (j);
                int idGlobal_j = pt_j->GetGlobalIndex ();
                LambdaOnPoint<Point> grad_phi_j = FEElement4Cell->GetGradPhi (j);
                LambdaOnPoint<real_t> phi_j = FEElement4Cell->GetPhi (j);


                for (ul_t edgeLocId = 0; edgeLocId < numOfEdgesOnCell; ++edgeLocId)
                {
                    Edge* edge = listEdgesOnCell->at (edgeLocId);
                    int idGlobal_edge = edge->GetGlobalIndex ();


                    if (tagedgevec->vec.at (ul_t (idGlobal_edge)) != static_cast<int>(tagToApply))
                        continue;

                    Edge* edgeOnCellRef = FEElement4Cell->GetEdge (edgeLocId);

                    Point* normal = edgeNormals->vec.at (ul_t (idGlobal_edge));
                    real_t prodscal = ((*edge->GetCentroid () - *cell->GetCentroid ()) | *normal);

                    if (prodscal < 0)
                        *normal = -1 * *normal;

                    FEBase* FEElement4Edge = festore.GetElementFor (edgeOnCellRef);
                    QuadStore::QuadObject* quad4Edge = quadstore.Get (FEElement4Edge);

                    real_t coeff1, coeff2;
                    Point v1, v2;
                    real_t value1, value2;
                    real_t pdetJac;
                    Matrix3x3_eig pJacInvT;

                    coeff1 = 0.;
                    coeff2 = 0.;

                    for (ul_t k = 0; k < quad4Edge->npts; ++k)
                    {
                        real_t w_k = quad4Edge->w [k];
                        //      Point pt_int_k = quad4Edge->pts [k];

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


void SolitonRoutines::ComputeDirichletNoSBM (Mesh* mesh, std::vector<Triplet_eig>* listTriplets, PlainVector_eig* F, real_t (*f)(Point, real_t), PHYS tagToApply, real_t hsz, real_t t)
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

    real_t max_aij = listTriplets->front ().value ();
    for (Triplet_eig triplet : *listTriplets)
        max_aij = std::max (max_aij, triplet.value ());

    real_t alpha = 20 / std::abs(hsz);

    INFOS << "Impose Dirichlet : tag " << to_string(tagToApply) << ", penal. coeff. = " << alpha << ENDLINE;
    int numCells = mesh->GetNumberOfCells ();

    auto tagedgevec = mesh->GetEdgesData ()->Get<int> (NAME_TAG_PHYSICAL);

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
        ul_t numOfEdgesOnCell = listEdgesOnCell->size ();
        std::vector<Point*>* listPointsOnCell = cell->GetPoints ();
        ul_t numOfPointsOnCell = listPointsOnCell->size ();

        FEBase* FEElement4Cell = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);

        for (ul_t i = 0; i < numOfPointsOnCell; ++i)
        {
            Point* pt_i = listPointsOnCell->at (i);
            int idGlobal_i = pt_i->GetGlobalIndex ();
            LambdaOnPoint<Point> grad_phi_i = FEElement4Cell->GetGradPhi (i);
            LambdaOnPoint<real_t> phi_i = FEElement4Cell->GetPhi (i);

            for (ul_t j = 0; j < numOfPointsOnCell; ++j)
            {
                Point* pt_j = listPointsOnCell->at (j);
                int idGlobal_j = pt_j->GetGlobalIndex ();
                LambdaOnPoint<Point> grad_phi_j = FEElement4Cell->GetGradPhi (j);
                LambdaOnPoint<real_t> phi_j = FEElement4Cell->GetPhi (j);


                for (ul_t edgeLocId = 0; edgeLocId < numOfEdgesOnCell; ++edgeLocId)
                {
                    Edge* edge = listEdgesOnCell->at (edgeLocId);
                    int idGlobal_edge = edge->GetGlobalIndex ();


                    if (tagedgevec->vec.at (ul_t (idGlobal_edge)) != static_cast<int>(tagToApply))
                        continue;

                    Edge* edgeOnCellRef = FEElement4Cell->GetEdge (edgeLocId);


                    FEBase* FEElement4Edge = festore.GetElementFor (edgeOnCellRef);
                    QuadStore::QuadObject* quad4Edge = quadstore.Get (FEElement4Edge);


                    real_t coeff1, coeff2, ck;
                    Point v1, v2;
                    Matrix3x3_eig pJacInvT;
                    real_t phi_i_eval, phi_j_eval;
                    Point grad_phi_i_eval, grad_phi_j_eval, normal_eval;
                    real_t u_D;

                    coeff1 = 0.;
                    coeff2 = 0.;
                    ck = 0.;
                    phi_i_eval = 0.;
                    phi_j_eval = 0.;
                    grad_phi_i_eval = {0., 0., 0.};
                    grad_phi_j_eval = {0., 0., 0.};
                    normal_eval = {0., 0., 0.};
                    u_D = 0.;

                    for (ul_t k = 0; k < quad4Edge->npts; ++k)
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
                        //      normal_eval = JacInvTEdge * JacInvTCell * *normal;
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



//void SolitonRoutines::ComputeDirichletSBM (Mesh* mesh, Mesh* object, std::vector<Triplet>* listTriplets, PlainVector* F, real_t (*f)(Point, real_t), real_t hsz, real_t t)
//{

// // * -<phi_j, grad(u).n>
// // * - <grad(phi_j).n, u - u_D>
// // * + <alpha * phi_j, u - u_D>

// // decouple
// // * A : -<phi_j, grad(u).n>
// // * A : - <grad(phi_j).n, u>
// // * A : + <alpha * phi_j, u>
// //
// // * B : - <grad(phi_j).n, u_D>
// // * B : + <alpha * phi_j, u_D>

// FEStore festore;
// QuadStore quadstore;

// real_t max_aij = listTriplets->front ().value ();
// for (Triplet_eig triplet : *listTriplets)
//  max_aij = std::max (max_aij, triplet.value ());

// real_t alpha = 20 / std::abs(hsz);

// INFOS << "Impose Dirichlet : name of vector of tag surrogate " << object->GetName ()+ NAME_TAG_INTERSECTION << ", penal. coeff. = " << alpha << ENDLINE;
// int numCells = mesh->GetNumberOfCells ();

// HetContainer<int>::Array* tagedgevec = mesh->GetEdgesData ()->Get<int> (object->GetName ()+ NAME_TAG_INTERSECTION);

// if (tagedgevec == nullptr)
// {
//  ERROR << "no array with the name : " << object->GetName ()+ NAME_TAG_INTERSECTION << "." << BLINKRETURN << ENDLINE;
//  return;
// }

// HetContainer<Point*>::Array* dvec = mesh->GetPointsData ()->Get<Point*> (object->GetName () + NAME_DISPLACEMENT_VECTOR);
// if (dvec == nullptr)
// {
//  ERROR << "no array with the name : " << object->GetName () + NAME_DISPLACEMENT_VECTOR << "." << BLINKRETURN << ENDLINE;
//  return;
// }

// FELocalInfos locEdge, locCell;


// for (int cellId = 0; cellId < numCells; ++cellId)
// {
//  Cell* cell = mesh->GetCell (cellId);

//  if (cell->GetTypeVTK () == VTK_LINE || cell->GetTypeVTK () == VTK_VERTEX)
//   continue;

//  std::vector<Edge*>* listEdgesOnCell = cell->GetEdges ();
//  ul_t numOfEdgesOnCell = listEdgesOnCell->size ();
//  std::vector<Point*>* listPointsOnCell = cell->GetPoints ();
//  ul_t numOfPointsOnCell = listPointsOnCell->size ();

//  FEBase* FEElement4Cell = festore.GetElementFor (cell, FE_CLASS_TYPE::LAGRANGE);

//  for (ul_t i = 0; i < numOfPointsOnCell; ++i)
//  {
//   Point* pt_i = listPointsOnCell->at (i);
//   int idGlobal_i = pt_i->GetGlobalIndex ();
//   LambdaOnPoint<Point> grad_phi_i = FEElement4Cell->gradPhi [i];
//   LambdaOnPoint<real_t> phi_i = FEElement4Cell->Phi [i];

//   for (ul_t j = 0; j < numOfPointsOnCell; ++j)
//   {
//    Point* pt_j = listPointsOnCell->at (j);
//    int idGlobal_j = pt_j->GetGlobalIndex ();
//    LambdaOnPoint<Point> grad_phi_j = FEElement4Cell->gradPhi [j];
//    LambdaOnPoint<real_t> phi_j = FEElement4Cell->Phi [j];


//    for (ul_t edgeLocId = 0; edgeLocId < numOfEdgesOnCell; ++edgeLocId)
//    {
//     Edge* edge = listEdgesOnCell->at (edgeLocId);
//     int idGlobal_edge = edge->GetGlobalIndex ();


//     if (tagedgevec->vec.at (ul_t (idGlobal_edge)) != tagToApply)
//      continue;

//     Edge* edgeOnCellRef = FEElement4Cell->edges [edgeLocId];

////     Point* normal = edgeNormals->vec.at (ul_t (idGlobal_edge));
////     real_t prodscal = ((*edge->GetCentroid () - *cell->GetCentroid ()) | *normal);

////     if (prodscal < 0)
////      *normal = -1 * *normal;

//     FEBase* FEElement4Edge = festore.GetElementFor (edgeOnCellRef);
//     QuadStore::QuadObject* quad4Edge = quadstore.Get (FEElement4Edge);


//     real_t coeff1, coeff2, ck;
//     Point v1, v2;
//     Matrix3x3 pJacInvT;
//     real_t phi_i_eval, phi_j_eval;
//     Point grad_phi_i_eval, grad_phi_j_eval, normal_eval;
//     real_t u_D;

//     coeff1 = 0.;
//     coeff2 = 0.;
//     ck = 0.;
//     phi_i_eval = 0.;
//     phi_j_eval = 0.;
//     grad_phi_i_eval = {0., 0., 0.};
//     grad_phi_j_eval = {0., 0., 0.};
//     normal_eval = {0., 0., 0.};
//     u_D = 0.;

//     for (ul_t k = 0; k < quad4Edge->npts; ++k)
//     {

//      Point pt_k_oncell = FEElement4Edge->TransformRefToEle (&quad4Edge->pts [k]);
//      Point pt_k_real = FEElement4Cell->TransformRefToEle (&pt_k_oncell);

//      FEElement4Edge->LocalCompute (&quad4Edge->pts [k], &locEdge);
//      FEElement4Cell->LocalCompute (&pt_k_oncell, &locCell);

//      ck = locEdge.detJac * locCell.detJac * quad4Edge->w [k];

//      phi_i_eval = phi_i (&pt_k_oncell);
//      phi_j_eval = phi_j (&pt_k_oncell);
//      grad_phi_i_eval = locEdge.JacInvT * locCell.JacInvT * grad_phi_i (&pt_k_oncell);
//      grad_phi_j_eval = locEdge.JacInvT * locCell.JacInvT * grad_phi_j (&pt_k_oncell);
//      //      normal_eval = JacInvTEdge * JacInvTCell * *normal;
//      normal_eval = {0, 1, 0};

//      u_D = f(pt_k_real, t);

//      // * A : -<phi_j, grad(phi_i).n>

//      coeff1 -= ck * phi_j_eval * (grad_phi_i_eval | normal_eval);

//      // * A : - <grad(phi_j).n, phi_i>

//      coeff1 -= ck * (grad_phi_j_eval | normal_eval) * phi_i_eval;

//      // * A : + <alpha * phi_j, phi_i>

//      coeff1 += ck * alpha * phi_j_eval * phi_i_eval;

//      if (i == 0)
//      {
//       // Second membre

//       // * - <grad(phi_j).n, u_D>
//       coeff2 -= ck * (grad_phi_j_eval | normal_eval) * u_D;

//       // * + <alpha * phi_j, u_D>
//       coeff2 += ck * alpha * phi_j_eval * u_D;
//      }

////      if (i == 0 && idGlobal_i == 283)
////       INFOS << "(i) 283 ck = " << coeff1 << " " << cellId << ENDLINE;
////      if (j == 0 && idGlobal_j == 283)
////       INFOS << "(j) 283 ck = " << coeff1 << " " << cellId << ENDLINE;

////      if (i == 0 && idGlobal_i == 279)
////       INFOS << "(i) 279 ck = " << coeff1 << " " << cellId << ENDLINE;
////      if (j == 0 && idGlobal_j == 279)
////       INFOS << "(j) 279 ck = " << coeff1 << " " << cellId << ENDLINE;
//     }

//     listTriplets->push_back({idGlobal_i, idGlobal_j, coeff1});

//     if (i == 0)
//      F->coeffRef (idGlobal_j) += coeff2;
//    }
//   }
//  }
// }
// return;
//}

void SolitonRoutines::PostProcess (MeshStorage *store, std::vector<SolitonFEContainer*> *containers, std::vector<PlainVector_eig*> *sols, std::function<real_t(Point, real_t)> fun)
{
    std::vector<INTER> list_tags_inter = {};
    for (SolitonFEContainer* cont : *containers)
        list_tags_inter.push_back (cont->GetInterTag ());

    // ul_t numPoints = static_cast<ul_t>(store->GetMainMesh ()->GetNumberOfPoints ());
    // ul_t numCells = static_cast<ul_t>(store->GetMainMesh ()->GetNumberOfCells ());
    ul_t numEdges = static_cast<ul_t>(store->GetMainMesh ()->GetNumberOfEdges ());

    HetContainer<int>::Array* tagSurEdge = store->GetMainMesh ()->GetEdgesData ()->Get<int> (NAME_TAG_INTERSECTION);

    for (ul_t idvec = 0; idvec < sols->size (); ++idvec)
    {
        STATUS << "postprocess on solver with tag " << COLOR_RED << to_string(list_tags_inter.at (idvec)) << ENDLINE;

        PlainVector_eig *solnum = sols->at (idvec);

        PlainVector_eig sol_ana;
        FunToVec (&sol_ana, store->GetMainMesh (), fun);

        if (list_tags_inter.at (idvec) != INTER::DEFAULT)
            for (ul_t idEdge = 0; idEdge < numEdges; ++idEdge)
            {
                if (tagSurEdge->vec.at (idEdge) == int(list_tags_inter.at (idvec)) ||
                        tagSurEdge->vec.at (idEdge) == int(INTER::MIXED))
                    continue;

                for (Point* pt : *store->GetMainMesh ()->GetEdge (int(idEdge))->GetPoints ())
                {
                    bool pass = false;
                    for (Cell* cell : pt->GetLinkedCell ())
                    {
                        if (cell->GetCat () == CAT_CELL_EDGE::CELL)
                            continue;

                        if (tagSurEdge->vec.at (static_cast<ul_t> (cell->GetGlobalIndex ())) == int(INTER::MIXED))
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

        if (tagsolver == "_" + to_string(INTER::DEFAULT))
            tagsolver = "";

        std::vector<real_t> solstd = PlainVector2Vector (solnum);
        store->GetMainMesh ()->GetPointsData ()->Add<real_t> ("sol_num" + tagsolver, solstd);

        solstd = PlainVector2Vector (&sol_ana);
        store->GetMainMesh ()->GetPointsData ()->Add<real_t> ("sol_ana" + tagsolver, solstd);

        PlainVector_eig erroabs = GetErrorAbs (store->GetMainMesh (), &sol_ana, solnum);
        solstd = PlainVector2Vector (&erroabs);
        store->GetMainMesh ()->GetPointsData ()->Add<real_t> ("err_abs" + tagsolver, solstd);

        //  PlainVector_eig errorelapercent = GetErrorRelaPercent (store->GetMainMesh (), &sol_ana, solnum);
        //  solstd = PlainVector2Vector (&errorelapercent);
        //  store->GetMainMesh ()->GetPointsData ()->Add<real_t> ("err_rela_percent" + tagsolver, solstd);

        GetErrorl1 (store->GetMainMesh (), &sol_ana, solnum);
        GetErrorl2 (store->GetMainMesh (), &sol_ana, solnum);
        GetErrorlinf (store->GetMainMesh (), &sol_ana, solnum);
        GetErrorRela (store->GetMainMesh (), &sol_ana, solnum);

        ENDFUN;
    }
}

