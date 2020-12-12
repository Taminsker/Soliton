#include <Soliton>

#include "functions.hpp"

const real_t beg_hsz    = 0.1;
const real_t end_hsz    = 99E-4;
const real_t step_hsz   = 1E-3;
const int    size_space = 20;

int
main (int argc, char * argv [])
{
    if (argc < 3)
    {
        ERROR << "Use : " << argv [0] << " input.slt out.dat" << ENDLINE;
        return EXIT_FAILURE;
    }

    bool enable_sbm        = false;
    bool enable_collect    = false;
    bool enable_view_all   = false;
    bool enable_view_once  = false;
    bool enable_synthetize = false;

    Eigen::setNbThreads (4);
    Eigen::initParallel ();

    SolitonRoutines::PrintExecutableName ("Poisson/Soliton");
    SolitonRoutines::PrintMacros ();

    STATUS << "enable collect   " << std::boolalpha << enable_collect << ENDLINE;
    STATUS << "enable sbm       " << std::boolalpha << enable_sbm << ENDLINE;
    STATUS << "enable synth     " << std::boolalpha << enable_synthetize << ENDLINE;
    STATUS << "enable view all  " << std::boolalpha << enable_view_all << ENDLINE;
    STATUS << "enable view once " << std::boolalpha << enable_view_once << ENDLINE << std::fixed;

    InputDatStruct inputdatfile;

    ParseInputDatFile (&inputdatfile, argv [1]);
    inputdatfile.grid_x_m = 0.0;
    inputdatfile.grid_x_p = 1.0;
    inputdatfile.grid_y_m = 0.0;
    inputdatfile.grid_y_p = 1.0;

    std::vector<real_t> err_l2, err_l1, err_linf, err_rela, hpres;

    COUT << std::setw (size_space) << "hsize" << std::setw (size_space) << "l1-error";
    COUT << std::setw (size_space) << "l2-error" << std::setw (size_space) << "linf-error";
    COUT << std::setw (size_space) << "rela-error" << std::setw (size_space) << "%" << ENDLINE;

    std::ofstream outfile (argv [2]);
    outfile << "# xmin " << inputdatfile.grid_x_m << std::endl;
    outfile << "# xmax " << inputdatfile.grid_x_p << std::endl;
    outfile << "# ymin " << inputdatfile.grid_y_m << std::endl;
    outfile << "# ymax " << inputdatfile.grid_y_p << std::endl;

    outfile << "#" << std::setw (size_space) << "hsize" << std::setw (size_space) << "l1-error";
    outfile << std::setw (size_space) << "l2-error" << std::setw (size_space) << "linf-error";
    outfile << std::setw (size_space) << "rela-error" << std::setw (size_space) << "%" << std::endl;

    for (real_t hsz = beg_hsz; hsz > end_hsz; hsz -= step_hsz)
    {
        hpres.push_back (hsz);
        inputdatfile.hsize = hsz;

        std::cout.setstate (std::ios_base::failbit);

        MeshStorage store;
        SolitonRoutines::Generation (&store, &inputdatfile);

        int                        numPoints    = store.GetMainMesh ()->GetNumberOfPoints ();
        int                        numEdges     = store.GetMainMesh ()->GetNumberOfEdges ();
        PointsData *               ptdata       = store.GetMainMesh ()->GetPointsData ();
        HetContainer<int>::Array * tagInterEdge = store.GetMainMesh ()->GetEdgesData ()->Get<int> (NAME_TAG_INTERSECTION);

        std::vector<SolitonFEContainer *> list_containers;
        std::vector<DenseVector>          identities;
        std::vector<std::string>          tags_str;

        real_t error_l1   = 0.;
        real_t error_l2   = 0.;
        real_t error_rela = 0.;
        real_t error_linf = 0.;

        if (enable_sbm)
        {
            list_containers.push_back (new SolitonFEContainer (&store, INTER::IN));
            list_containers.push_back (new SolitonFEContainer (&store, INTER::OUT));

            identities.push_back (DenseVector::Zero (numPoints));  // IN
            identities.push_back (DenseVector::Zero (numPoints));  // OUT

            tags_str.push_back (to_string (INTER::IN));
            tags_str.push_back (to_string (INTER::OUT));
        }
        else
        {
            list_containers.push_back (new SolitonFEContainer (&store, INTER::DEFAULT));

            identities.push_back (DenseVector::Ones (numPoints));

            tags_str.push_back ("");
        }

        for (ul_t id_inter = 0; id_inter < list_containers.size (); ++id_inter)
        {
            SolitonFEContainer * container = list_containers.at (id_inter);

            // GRAD_GRAD
            SolitonFEItem<ITEM_T::GRAD_GRAD> ggI (container, PHYS::DOMAIN);
            ggI.SetFunction2Apply (&fun_dirichlet);
            container->AddItem (ggI);

            // SECOND_MEMBER
            SolitonFEItem<ITEM_T::SECOND_MEMBER> scI (container, PHYS::DOMAIN);
            scI.SetFunction2Apply (&fun_secondmember);
            container->AddItem (scI);

            // DIRICHLET_BOUNDARY : INLET
            SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> dIi (container, PHYS::INLET);
            dIi.SetFunction2Apply (&fun_dirichlet);
            container->AddItem (dIi);

            // DIRICHLET_BOUNDARY : OUTLET
            SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> dIo (container, PHYS::OUTLET);
            dIo.SetFunction2Apply (&fun_dirichlet);
            container->AddItem (dIo);

            // DIRICHLET_BOUNDARY : WALL
            SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY> dIw (container, PHYS::WALL);
            dIw.SetFunction2Apply (&fun_dirichlet);
            container->AddItem (dIw);

            if (enable_sbm)
            {
                // DIRICHLET_BOUNDARY_SBM
                SolitonFEItem<ITEM_T::DIRICHLET_BOUNDARY_SBM> dIsbm (container, PHYS::DEFAULT);
                dIsbm.SetFunction2Apply (&fun_dirichlet);
                container->AddItem (dIsbm);
            }

            container->SetCoeffPen (inputdatfile.penal);
            container->SetPowPen (inputdatfile.powpenalty);
            container->SetCollectContribution (enable_collect);

            // Identy part
            for (int edgeId = 0; edgeId < numEdges; ++edgeId)
            {
                INTER tag_edge = static_cast<INTER> (tagInterEdge->vec.at (static_cast<ul_t> (edgeId)));

                for (Point * pt : *store.GetMainMesh ()->GetEdge (edgeId)->GetPoints ())
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

        for (SolitonFEContainer * container : list_containers)
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
       * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER
       * =>
       * MAT_0*U = SECMEMBER
       */

        for (SolitonFEContainer * container : list_containers)
        {
            if (!enable_synthetize)
                BEGIN << "Compute on INTER::" << to_string (container->GetInterTag ()) << ENDLINE;

            std::vector<SparseMatrix> listMats;
            DenseVector               sec;
            container->ComputeSecMemberAndMats (&listMats, &sec, 0x0);

            SparseMatrix matrix (numPoints, numPoints);
            DenseVector  secMember = DenseVector::Zero (numPoints);

            matrix    = listMats.at (0);
            secMember = sec;

            DenseVector solnum;
            Solver::AutoDeduceBest (&matrix, &secMember, &solnum, !enable_synthetize);

            container->GetQueue (0)->New (solnum);
        }

        // Post-process

        std::vector<real_t> solstd;
        DenseVector     sol_ana_syn, sol_num_syn, err_abs_syn;
        sol_num_syn.setZero (numPoints);
        err_abs_syn.setZero (numPoints);

        FunToVec (&sol_ana_syn, store.GetMainMesh (), &fun_analytic);
        solstd = PlainVector2Vector (&sol_ana_syn);
        ptdata->Add<real_t> ("sol_ana", solstd);

        for (ul_t id_inter = 0; id_inter < list_containers.size (); ++id_inter)
        {
            SolitonFEContainer * container = list_containers.at (id_inter);

            if (!enable_synthetize)
                BEGIN << "Postprocess on INTER::" << to_string (container->GetInterTag ()) << ENDLINE;

            std::string       tag_str_inter = to_string (container->GetInterTag ());
            DenseVector * sol_num       = container->GetQueue (0)->GetNumber (0);

            *sol_num = sol_num->cwiseProduct (identities.at (id_inter));
            sol_num_syn += *sol_num;
            DenseVector sol_ana = sol_ana_syn.cwiseProduct (identities.at (id_inter));

            solstd = PlainVector2Vector (sol_num);

            if (!enable_synthetize)
                ptdata->Add<real_t> ("sol_num" + tag_str_inter, solstd);

            DenseVector erroabs = GetErrorAbs (store.GetMainMesh (), &sol_ana, sol_num, !enable_synthetize);
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

        error_l1   = GetErrorl1 (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);
        error_l2   = GetErrorl2 (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);
        error_linf = GetErrorlinf (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);
        error_rela = GetErrorRela (store.GetMainMesh (), &sol_ana_syn, &sol_num_syn);

        err_l1.push_back (error_l1);
        err_l2.push_back (error_l2);
        err_linf.push_back (error_linf);
        err_rela.push_back (error_rela);

        ENDFUN;

        std::cout.clear ();

        COUT << std::scientific << FLUSHLINE;
        COUT << COLOR_YELLOW << std::setw (size_space) << hsz << COLOR_RED;
        COUT << std::setw (size_space) << err_l1.back ();
        COUT << std::setw (size_space) << err_l2.back ();
        COUT << std::setw (size_space) << err_linf.back ();
        COUT << std::setw (size_space) << err_rela.back ();
        COUT << std::fixed << FLUSHLINE;
        COUT << std::setw (size_space) << err_rela.back () * 100. << "%" << ENDLINE;

        outfile << std::scientific << std::flush;
        outfile << std::setw (size_space) << hsz;
        outfile << std::setw (size_space) << err_l1.back ();
        outfile << std::setw (size_space) << err_l2.back ();
        outfile << std::setw (size_space) << err_linf.back ();
        outfile << std::setw (size_space) << err_rela.back ();
        outfile << std::fixed << std::endl;
    }

    COUT << std::scientific << FLUSHLINE;
    std::cout << "err_l1 = [ ";
    for (ul_t i = 0; i < err_l1.size () - 1; ++i)
        std::cout << err_l1.at (i) << ", ";
    std::cout << err_l1.back () << " ];" << std::endl;

    std::cout << "err_l2 = [ ";
    for (ul_t i = 0; i < err_l2.size () - 1; ++i)
        std::cout << err_l2.at (i) << ", ";
    std::cout << err_l2.back () << " ];" << std::endl;

    std::cout << "err_linf = [ ";
    for (ul_t i = 0; i < err_linf.size () - 1; ++i)
        std::cout << err_linf.at (i) << ", ";
    std::cout << err_linf.back () << " ];" << std::endl;

    std::cout << "err_rela = [ ";
    for (ul_t i = 0; i < err_rela.size () - 1; ++i)
        std::cout << err_rela.at (i) << ", ";
    std::cout << err_rela.back () << " ];" << std::endl;

    std::cout << "h = [ ";
    for (ul_t i = 0; i < hpres.size () - 1; ++i)
        std::cout << hpres.at (i) << ", ";
    std::cout << hpres.back () << " ];" << std::endl;

    return EXIT_SUCCESS;
}
