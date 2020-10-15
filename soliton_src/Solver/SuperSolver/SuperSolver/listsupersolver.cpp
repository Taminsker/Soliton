
#include <Solver/SuperSolver/SuperSolver/supersolver.h>
#include <Solver/Method/method.h>

template<>
void SuperSolver::InternalCompute<SCH_T::NO_TIME> ()
{
    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     *  =>
     * MAT_0*U = SECMEMBER*U
     */

    this->Print ();

    std::vector<SparseMatrix> listMats;
    PlainVector sec;

    this->ComputeSecMemberAndMats (&listMats, &sec);

    int numPoints = m_store->mesh->GetNumberOfPoints ();
    SparseMatrix matrix(numPoints, numPoints);
    PlainVector secMember = PlainVector::Zero (numPoints);

    matrix = listMats.at (0);
    //    ERROR << "LISTITEM" << ENDLINE;
    //    std::cout << matrix << std::endl;
    secMember = sec;
    //        INFOS << sec.transpose ()<< ENDLINE;

    /*******************************/

    //    auto matTemp = *m_items [0][0]->CastSparseMatrix ();
    //    ERROR << SEPARATOR << "MATRIX LISTITEM MINUS GRAD" << SEPARATOR << ENDLINE;
    //    std::cout << Eigen::MatrixXd(matrix - *m_items [0][0]->CastSparseMatrix ()) << ENDLINE;
    //    ERROR << SEPARATOR << "NEW" << SEPARATOR << ENDLINE;

    //    matrix = *m_items [0][0]->CastSparseMatrix ();
    //    secMember.setZero ();
    //    secMember = *m_items [0][1]->CastPlainVector ();

    ////    INFOS << *m_items [0][1] << ENDLINE;

    //    auto g = m_items.at (0).at (2)->m_fun;

    //    for(int i = 0; i < numPoints; ++i)
    //    {
    //        PHYS tag = PHYS(m_store->mesh->GetPointsData ()->GetIntArrays ()->Get (NAME_PHYS)->vec.at (std::size_t (i)));

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

    //    ERROR << "MATRIX AFTER OLD DIRICHLET" << ENDLINE;
    //    std::ofstream outfile("txt_matrix.txt");
    //    outfile << matrix << std::endl;
    //    outfile.close ();

    //    ERROR << secMember.transpose () << ENDLINE;


    //    /******************************/

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    return m_queue->New (solution);
}


template<>
[[deprecated("Not set...")]]
void
SuperSolver::InternalCompute<SCH_T::EULER_EXPLICIT> ()
{
    this->Print ();

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
                sec += coeffs.at (id) * listMats.at (idDer) * *m_queue->GetNumber (idx.at (id));
    }

    matrix += listMats.at (0);
    secMember += sec;

    PlainVector solution;
    Solver::AutoDeduceBest (&matrix, &secMember, &solution);

    return m_queue->New (solution);
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

    this->Print ();

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

    return m_queue->New (solution);
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

    this->Print ();

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

    return m_queue->New (solution);
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

    this->Print ();

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

    return m_queue->New (solution);
}

template<>
[[deprecated("Not set...")]]
void
SuperSolver::InternalCompute<SCH_T::NEWMARK_SCHEME> ()
{
    /**
     * Dt_2*MAT_2*U + Dt_1*MAT_1*U + MAT_0*U = SECMEMBER*U
     */

    /**
      * Newmark-beta method (based on df temporal)
      * u'_{n+1} = u'_{n} + dt * u"_{alph}
      **** u"_{alpha} = (1-alpha) * u"_{n} + alpha * u"_{n}
      * u_{n+1} = u_{n} + dt * u'_{n} + 0.5 * dt^2 * u"_{beta}
      **** u"_{beta} = (1-2*beta) * u"_{n} + 2 * beta * u"_{n+1}
      **/

    this->Print ();

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

    return m_queue->New (solution);
}
