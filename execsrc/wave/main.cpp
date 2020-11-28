#include <Soliton>
#include "functions.h"

int main(int argc, char *argv[])
{
    coeffA          = 0.1;
    coeffT          = 0.5;
    coeffLambda     = 0.2;

    real_t celerity = coeffLambda / coeffT;

    real_t dt = coeffT / 10.;

    if (argc < 3)
    {
        ERROR << "Use : " << argv[0] << " input.slt nIter [options...]" << ENDLINE;
        ERROR << "options : " << ENDLINE;
        ERROR << "  -collect            : collect items contributions" << ENDLINE;
        ERROR << "  -sbm                : sbm on objects interfaces" << ENDLINE;
        ERROR << "  -synth              : synthetize results" << ENDLINE;
        ERROR << "  -view-all           : view all compute steps" << ENDLINE;
        ERROR << "  -view-once          : view the conf step" << ENDLINE;
        return EXIT_FAILURE;
    }

    int maxIter = std::stoi(argv[2]);

    bool enable_sbm = false;
    bool enable_collect = false;
    bool enable_view_all = false;
    bool enable_view_once = false;
    bool enable_synthetize = false;

    for (int i = 2; i < argc; ++i)
    {
        enable_collect      = (std::string(argv[i]) == "-collect")      || enable_collect;
        enable_sbm          = (std::string(argv[i]) == "-sbm")          || enable_sbm;
        enable_synthetize   = (std::string(argv[i]) == "-synth")        || enable_synthetize;
        enable_view_all     = (std::string(argv[i]) == "-view-all")     || enable_view_all;
        enable_view_once    = (std::string(argv[i]) == "-view-once")    || enable_view_once;
    }

    Eigen::setNbThreads(4);
    Eigen::initParallel();

    SolitonRoutines::PrintExecutableName ("Poisson/Soliton");
    SolitonRoutines::PrintMacros ();

    STATUS << "enable collect   " << std::boolalpha << enable_collect << ENDLINE;
    STATUS << "enable sbm       " << std::boolalpha << enable_sbm << ENDLINE;
    STATUS << "enable synth     " << std::boolalpha << enable_synthetize << ENDLINE;
    STATUS << "enable view all  " << std::boolalpha << enable_view_all << ENDLINE;
    STATUS << "enable view once " << std::boolalpha << enable_view_once << ENDLINE << std::fixed;
    ENDFUN;

    InputDatStruct inputdatfile;
    ParseInputDatFile (&inputdatfile, argv [1]);

    MeshStorage store;
    SolitonRoutines::Generation (&store, &inputdatfile);

    int numPoints                           = store.GetMainMesh ()->GetNumberOfPoints ();
    int numEdges                            = store.GetMainMesh ()->GetNumberOfEdges ();
    PointsData *ptdata                      = store.GetMainMesh ()->GetPointsData ();
    HetContainer<int>::Array* tagInterEdge  = store.GetMainMesh ()->GetEdgesData ()->Get<int>(NAME_TAG_INTERSECTION);

    std::vector<SolitonFEContainer*>  list_containers;
    std::vector<PlainVector_eig>    identities;
    std::vector<std::string>        tags_str;

    if (enable_sbm)
    {
        list_containers.push_back (new SolitonFEContainer (&store, INTER::IN));
        list_containers.push_back (new SolitonFEContainer (&store, INTER::OUT));

        identities.push_back (PlainVector_eig::Zero (numPoints)); // IN
        identities.push_back (PlainVector_eig::Zero (numPoints)); // OUT

        tags_str.push_back (to_string(INTER::IN));
        tags_str.push_back (to_string(INTER::OUT));
    }
    else
    {
        list_containers.push_back (new SolitonFEContainer (&store, INTER::DEFAULT));

        identities.push_back (PlainVector_eig::Ones (numPoints));

        tags_str.push_back ("");
    }


    for (ul_t id_inter = 0; id_inter < list_containers.size (); ++id_inter)
    {
        SolitonFEContainer *container = list_containers.at (id_inter);

//        // PHI_PHI
//        SolitonFEItem<ITEM_T::PHI_PHI> ppI (container, PHYS::DOMAIN);
//        ppI.SetTemporalOrderDerivative (2);
//        container->AddItem(ppI);

        // GRAD_GRAD
        SolitonFEItem<ITEM_T::GRAD_GRAD> ggI (container, PHYS::DOMAIN);
        container->AddItem(ggI);

        // SECOND_MEMBER
        SolitonFEItem<ITEM_T::SECOND_MEMBER> scI (container, PHYS::DOMAIN);
        scI.SetFunction2Apply (&fun_secondmember);
        container->AddItem(scI);

        // DIRICHLET_BOUNDARY : INLET
        SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> dIi (container, PHYS::INLET);
        dIi.SetFunction2Apply (&fun_dirichlet);
        dIi.SetVariableOverTimeOn ();
        container->AddItem(dIi);

        // DIRICHLET_BOUNDARY : OUTLET
        SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> dIo (container, PHYS::OUTLET);
        dIo.SetFunction2Apply (&fun_dirichlet);
        dIo.SetVariableOverTimeOn ();
        container->AddItem(dIo);

        //        // NEUMANN_BOUNDARY : WALL
        //        SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> nIw (container, PHYS::WALL);
        //        nIw.SetFunction2Apply (&fun_neumann);
        //        container->AddItem(nIw);

        // DIRICHLET_BOUNDARY : WALL
        SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> dIw (container, PHYS::WALL);
        dIw.SetFunction2Apply (&fun_dirichlet);
        dIw.SetVariableOverTimeOn ();
        container->AddItem(dIw);

        if (enable_sbm)
        {
            // DIRICHLET_BOUNDARY_SBM
            SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY_SBM> dIsbm (container, PHYS::DEFAULT);
            dIsbm.SetFunction2Apply (&fun_dirichlet);
            dIsbm.SetVariableOverTimeOn ();
            container->AddItem(dIsbm);
        }

        container->SetCoeffPen (inputdatfile.penal);
        container->SetPowPen (inputdatfile.powpenalty);
        container->SetCollectContribution (enable_collect);

        // Identy part
        for (int edgeId = 0; edgeId < numEdges; ++edgeId)
        {
            INTER tag_edge = static_cast<INTER>(tagInterEdge->vec.at (static_cast<ul_t>(edgeId)));

            for (Point* pt : *store.GetMainMesh ()->GetEdge (edgeId)->GetPoints ())
            {
                if (container->GetInterTag () == INTER::DEFAULT)
                    identities.at (id_inter).coeffRef (pt->GetGlobalIndex ()) = 1.;
                else if (tag_edge == container->GetInterTag ())
                    identities.at (id_inter).coeffRef (pt->GetGlobalIndex ()) = 1.;
                else if (tag_edge == INTER::MIXED)
                    identities.at (id_inter).coeffRef (pt->GetGlobalIndex ()) = 0.5;
            }
        }
    }


    for (SolitonFEContainer* container : list_containers)
    {
        container->SetView (enable_view_once);
        container->SetForceComputeOn ();

        container->Configure ();

        container->SetView (enable_view_all);
        container->SetForceComputeOff ();
        container->SetCollectContributionOff ();
    }

    // Compute part (no_time)

    /**
         * Dt_2*Mat_2*U + Dt_1*Mat_1*U + Mat_0*U = SECMEMBER
         */

    /**
          * Newmark-beta method (based on df temporal)
          *
          * u'_{n+1} = u'_{n} + dt * u"_{alpha}
          **** u"_{alpha} = (1-alpha) * u"_{n+1} + alpha * u"_{n}
          *
          * u_{n+1} = u_{n} + dt * u'_{n} + 0.5 * dt^2 * u"_{beta}
          **** u"_{beta} = (1-2*beta) * u"_{n} + 2 * beta * u"_{n+1}
          **/

    real_t alpha = 0.5;
    real_t beta = 0.0;
    INFOS << "Newmark scheme with alpha = " << alpha << " and beta = " << beta << ENDLINE;

    // init
    for (SolitonFEContainer* container : list_containers)
    {
        PlainVector_eig *u_0_0 = container->GetQueue (0)->GetNumber (0);
        PlainVector_eig *u_1_0 = container->GetQueue (1)->GetNumber (0);
        PlainVector_eig *u_2_0 = container->GetQueue (2)->GetNumber (0);

        FunToVec (u_0_0, store.GetMainMesh (), &fun_analytic, 0x0);
        FunToVec (u_1_0, store.GetMainMesh (), &fun_analytic_1, 0x0);
        FunToVec (u_2_0, store.GetMainMesh (), &fun_analytic_2, 0x0);

        std::vector<real_t> solstd;

        solstd = PlainVector2Vector (u_0_0);
        ptdata->Add<real_t> ("sol_ana", solstd);
        ptdata->Add<real_t> ("sol_num", solstd);

        solstd = std::vector<real_t>(static_cast<ul_t>(numPoints), 0x0);
        ptdata->Add<real_t> ("err_abs", solstd);


        // Write part
        WriteVTKWithCells (store.GetMainMesh (), std::to_string (0), false);
        WriteVTKWithEdges (store.GetMainMesh (), std::to_string (0), false);

        for (Mesh* object : *store.GetListOfObjects ())
        {
            WriteVTKWithCells (object, "", false);
            WriteVTKWithEdges (object, "", false);
        }
        ptdata->Remove<real_t> ("sol_num");
        ptdata->Remove<real_t> ("sol_ana");
        ptdata->Remove<real_t> ("err_abs");
    }

    for (int itTime = 1; itTime < maxIter; ++itTime)
    {
        real_t time = static_cast<real_t> (itTime) * dt;

        for (SolitonFEContainer* container : list_containers)
        {
            if (!enable_synthetize)
                BEGIN << "Compute on INTER::" << to_string(container->GetInterTag ()) << ENDLINE;

            std::vector<SparseMatrix_eig> listMats;
            PlainVector_eig sec;
            container->ComputeSecMemberAndMats (&listMats, &sec, time);

            // Init
            SparseMatrix_eig matrix(numPoints, numPoints);
            PlainVector_eig secMember = PlainVector_eig::Zero (numPoints);

            PlainVector_eig *u_0_n = container->GetQueue (0)->GetNumber (0);
            PlainVector_eig *u_1_n = container->GetQueue (1)->GetNumber (0);
            PlainVector_eig *u_2_n = container->GetQueue (2)->GetNumber (0);

            listMats.at (0) = celerity * celerity * listMats.at (0);
            listMats.at (2).setIdentity ();

            // ---------------------------------------------- //
            // Solve u_2_np1
            //
            // Mat_2*U_2_np1 = SecMember - Mat_1 * U_1_n - Mat_0 * U_0_n

//            matrix = listMats.at (2);
//            std::ofstream mat("matrix.txt");
//            mat << matrix;
//            mat.close ();

            secMember = sec - listMats.at (1) * *u_1_n - listMats.at (0) * *u_0_n;

            PlainVector_eig u_2_np1;
            Solver::AutoDeduceBest (&matrix, &secMember, &u_2_np1);
            container->GetQueue (2)->New (u_2_np1);

            // ---------------------------------------------- //
            // Solve u_1_np1
            //
            // U_1_np1 = U_1_n + dt * U_2_alpha;
            // with
            // U_2_alpha = alpha * U_2_n + (1. - alpha) * U_2_np1;

            PlainVector_eig u_2_alpha = alpha * *u_2_n + (1. - alpha) * u_2_np1;
            PlainVector_eig u_1_np1 = *u_1_n + dt * u_2_alpha;
            container->GetQueue (1)->New (u_1_np1);

            // ---------------------------------------------- //
            // Solve u_0_np1
            //
            // U_1_np1 = U_0_n + dt * U_1_n + 0.5 * dt*dt * U_2_beta;
            // with
            // U_2_beta = 2. * beta * U_2_n + (1. - 2. * beta) * U_2_np1;

            PlainVector_eig u_2_beta = 2. * beta * *u_2_n + (1. - 2. * beta) * u_2_np1;
            PlainVector_eig u_0_np1 = *u_0_n + dt * *u_1_n + 0.5 * dt * dt * u_2_beta;

            container->GetQueue (0)->New (u_0_np1);
        }

        // Post-process

        std::vector<real_t> solstd;
        PlainVector_eig sol_ana_syn, sol_num_syn, err_abs_syn;
        sol_num_syn.setZero (numPoints);
        err_abs_syn.setZero (numPoints);

        FunToVec (&sol_ana_syn, store.GetMainMesh (), &fun_analytic, time);
        solstd = PlainVector2Vector (&sol_ana_syn);
        ptdata->Add<real_t> ("sol_ana", solstd);


        for (ul_t id_inter = 0; id_inter < list_containers.size (); ++id_inter)
        {
            SolitonFEContainer *container = list_containers.at (id_inter);

            if (!enable_synthetize)
                BEGIN << "Postprocess on INTER::" << to_string(container->GetInterTag ()) << ENDLINE;

            std::string tag_str_inter = to_string(container->GetInterTag ());
            PlainVector_eig* sol_num = container->GetQueue (0)->GetNumber (0);

            *sol_num = sol_num->cwiseProduct(identities.at (id_inter));
            sol_num_syn += *sol_num;
            PlainVector_eig sol_ana = sol_ana_syn.cwiseProduct(identities.at (id_inter));

            solstd = PlainVector2Vector (sol_num);

            if (!enable_synthetize)
                ptdata->Add<real_t> ("sol_num" + tag_str_inter, solstd);

            PlainVector_eig erroabs = GetErrorAbs (store.GetMainMesh (), &sol_ana, sol_num, !enable_synthetize);
            err_abs_syn += erroabs.cwiseProduct (identities.at (id_inter));
            solstd = PlainVector2Vector (&erroabs);

            if (!enable_synthetize)
                ptdata->Add<real_t> ("err_abs" + tag_str_inter, solstd);

            GetErrorl1 (store.GetMainMesh (), &sol_ana, sol_num, !enable_synthetize);
            GetErrorl2 (store.GetMainMesh (), &sol_ana, sol_num, !enable_synthetize);
            GetErrorlinf (store.GetMainMesh (), &sol_ana, sol_num, !enable_synthetize);
            GetErrorRela (store.GetMainMesh (), &sol_ana, sol_num, !enable_synthetize);


            if (!enable_synthetize)
                ENDFUN;
        }

        BEGIN << "Postprocess global at time " << time << ENDLINE;
        GetErrorl1 (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);
        GetErrorl2 (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);
        GetErrorlinf (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);
        GetErrorRela (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);

        solstd = PlainVector2Vector (&sol_num_syn);
        ptdata->Add<real_t> ("sol_num", solstd);
        solstd = PlainVector2Vector (&err_abs_syn);
        ptdata->Add<real_t> ("err_abs", solstd);

        ENDFUN;


        // Write part
        WriteVTKWithCells (store.GetMainMesh (), std::to_string (itTime), false);
        WriteVTKWithEdges (store.GetMainMesh (), std::to_string (itTime), false);

        for (Mesh* object : *store.GetListOfObjects ())
        {
            WriteVTKWithCells (object, "", false);
            WriteVTKWithEdges (object, "", false);
        }
        ptdata->Remove<real_t> ("sol_num");
        ptdata->Remove<real_t> ("sol_ana");
        ptdata->Remove<real_t> ("err_abs");


    }

    return EXIT_SUCCESS;
}
