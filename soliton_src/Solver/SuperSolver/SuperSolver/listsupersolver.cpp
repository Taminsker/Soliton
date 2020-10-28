
#include <Solver/SuperSolver/SuperSolver/supersolver.h>
#include <Solver/SuperSolver/SuperItem/superitem.h>

#include <Solver/Method/method.h>

template<>
void SuperSolver::InternalCompute<SCH_T::NO_TIME> ()
{
    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     *  =>
     * MAT_0*U = SECMEMBER*U
     */

    std::vector<SparseMatrix> listMats;
    PlainVector sec;

    this->ComputeSecMemberAndMats (&listMats, &sec);

    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    matrix = listMats.at (0);
    secMember = sec;

//    /******************************/

//    auto matTemp = *m_list_items->at (0).at (0)->GetSparseMatrix ();
//    //        ERROR << SEPARATOR << "MATRIX LISTITEM MINUS GRAD" << SEPARATOR << ENDLINE;
//    //        std::cout << Eigen::MatrixXd(matrix - *m_list_items->at (0).at (0)->GetSparseMatrix ()) << ENDLINE;
//    //        ERROR << SEPARATOR << "NEW" << SEPARATOR << ENDLINE;

//    matrix = *m_list_items->at (0).at (0)->GetSparseMatrix ();
//    secMember.setZero ();
//    secMember = *m_list_items->at (0).at (1)->GetPlainVector ();

//    //    INFOS << *m_items [0][1] << ENDLINE;

//    auto g = m_list_items->at (0).at (0)->GetFun ();

//    for(int i = 0; i < numPoints; ++i)
//    {
//        PHYS tag = PHYS(m_store->mesh->GetPointsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL)->vec.at (std::size_t (i)));

//        if (tag == PHYS::WALL || tag == PHYS::INLET || tag == PHYS::OUTLET)
//        {
//            matrix.row (i) *= 0.;

//            auto temp = matrix.transpose (); // Pour pouvoir manipuler les colonnes de A
//            matrix = temp;

//            secMember -= g(*m_store->mesh->GetPoint (i), m_time) * matrix.row (i).transpose ();
//            matrix.row (i) *= 0.;
//            matrix.coeffRef (i,i) = 1.;
//            secMember.coeffRef (i) = g(*m_store->mesh->GetPoint (i), m_time);

//            temp = matrix.transpose ();
//            temp.pruned ();
//            matrix = temp;
//        }
//    }

//    //        ERROR << "MATRIX AFTER OLD DIRICHLET" << ENDLINE;
//    //        std::ofstream outfile("txt_matrix.txt");
//    //        outfile << matrix << std::endl;
//    //        outfile.close ();

//    //        ERROR << secMember.transpose () << ENDLINE;


//    /*****************************/

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    return m_queues.at (0)->New (solution);
}


template<>
[[deprecated("Not set...")]]
void
SuperSolver::InternalCompute<SCH_T::EULER_EXPLICIT> ()
{
    std::vector<SparseMatrix> listMats;
    PlainVector sec;

    this->ComputeSecMemberAndMats (&listMats, &sec);

    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     */

    std::vector<DF_ORDER_NAMES> orders = {DF_ORDER_NAMES::ORDER_1_CENTRAL, DF_ORDER_NAMES::ORDER_3_BACKWARD, DF_ORDER_NAMES::ORDER_2_BACKWARD};
    std::vector<int> idx;
    std::vector<double> coeffs;

    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    for (std::size_t idDer = 1; idDer < 3; ++idDer)
    {
        m_dfstore->Get (static_cast<int> (idDer), orders.at (idDer), &idx, &coeffs);
        for (std::size_t id = 0; id < idx.size (); ++id)
            if (idx.at (id) == 0)
                matrix += coeffs.at (id) * listMats.at (idDer);
            else
                sec += coeffs.at (id) * listMats.at (idDer) * *m_queues.at (0)->GetNumber (idx.at (id));
    }

    matrix += listMats.at (0);
    secMember += sec;

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    return m_queues.at (0)->New (solution);
}

template<>
[[deprecated("Not set...")]]
void
SuperSolver::InternalCompute<SCH_T::EULER_IMPLICIT> ()
{
    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     *  =>
     * MAT_0*U = SECMEMBER*U
     */

    std::vector<SparseMatrix> listMats;
    PlainVector sec;

    this->ComputeSecMemberAndMats (&listMats, &sec);

    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    matrix += listMats.at (0);
    secMember += sec;

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    return m_queues.at (0)->New (solution);
}

template<>
[[deprecated("Not set...")]]
void
SuperSolver::InternalCompute<SCH_T::TVD_RK_2> ()
{
    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     *  =>
     * MAT_0*U = SECMEMBER*U
     */

    std::vector<SparseMatrix> listMats;
    PlainVector sec;

    this->ComputeSecMemberAndMats (&listMats, &sec);

    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    matrix += listMats.at (0);
    secMember += sec;

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    ENDFUN;
    return m_queues.at (0)->New (solution);
}


template<>
[[deprecated("Not set...")]]
void
SuperSolver::InternalCompute<SCH_T::TVD_RK_4> ()
{
    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     *  =>
     * MAT_0*U = SECMEMBER*U
     */

    std::vector<SparseMatrix> listMats;
    PlainVector sec;

    this->ComputeSecMemberAndMats (&listMats, &sec);

    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    matrix += listMats.at (0);
    secMember += sec;

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    return m_queues.at (0)->New (solution);
}

template<>
void
SuperSolver::InternalCompute<SCH_T::NEWMARK_SCHEME> ()
{
    /**
     * Dt_2*Mat_2*U + Dt_1*Mat_1*U + Mat_0*U = SECMEMBER
     */

    /**
      * Newmark-beta method (based on df temporal)
      *
      * u'_{n+1} = u'_{n} + dt * u"_{alpha}
      **** u"_{alpha} = (1-alpha) * u"_{n} + alpha * u"_{n}
      *
      * u_{n+1} = u_{n} + dt * u'_{n} + 0.5 * dt^2 * u"_{beta}
      **** u"_{beta} = (1-2*beta) * u"_{n} + 2 * beta * u"_{n+1}
      **/

    double alpha = 0.5;
    double beta = 0.5;
    INFOS << "Newmark scheme with alpha = " << alpha << " and beta = " << beta << ENDLINE;

    // Build matrix and second member
    std::vector<SparseMatrix> listMats;
    PlainVector sec;
    this->ComputeSecMemberAndMats (&listMats, &sec);

    // Init
    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    PlainVector *u_0_n = m_queues.at (0)->GetNumber (0);
    PlainVector *u_1_n = m_queues.at (1)->GetNumber (0);
    PlainVector *u_2_n = m_queues.at (2)->GetNumber (0);

    // ---------------------------------------------- //
    // Solve u_2_np1
    //
    // Mat_2*U_2_np1 = SecMember - Mat_1 * U_1_n - Mat_0 * U_0_n

    matrix = listMats.at (2);
    secMember = sec - listMats.at (1) * *u_1_n - listMats.at (0) * *u_0_n;

    PlainVector u_2_np1;
    Solver::AutoDeduceBest (&matrix, &secMember, &u_2_np1);
    m_queues.at (2)->New (u_2_np1);

    // ---------------------------------------------- //
    // Solve u_1_np1
    //
    // U_1_np1 = U_1_n + dt * U_2_alpha;
    // with
    // U_2_alpha = alpha * U_2_n * (1. - alpha) * U_2_np1;

    PlainVector u_2_alpha = alpha * *u_2_n * (1. - alpha) * u_2_np1;
    PlainVector u_1_np1 = *u_1_n + m_dt * u_2_alpha;
    m_queues.at (1)->New (u_1_np1);

    // ---------------------------------------------- //
    // Solve u_0_np1
    //
    // U_1_np1 = U_0_n + dt * U_1_n + 0.5 * dt*dt * U_2_beta;
    // with
    // U_2_beta = 2. * beta * U_2_n * (1. - 2. * beta) * U_2_np1;

    PlainVector u_2_beta = 2. * beta * *u_2_n * (1. - 2. * beta) * u_2_np1;
    PlainVector u_0_np1 = *u_0_n + m_dt * *u_1_n + 0.5 * m_dt * m_dt * u_2_beta;

    return m_queues.at (0)->New (u_0_np1);
}
