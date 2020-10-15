//#include "upgradablesolverbase.h"
//#include <IO/parsers.h>

//UpgradableSolverBase::UpgradableSolverBase (Sto4Sol* store) :
//    m_store (store),
//    m_time (0.),
//    m_coeff_pen (1.),
//    m_dt (1.),
//    m_collectcontributions (false),
//    m_items (unstd::matrix<ItemSolverBase*> (UPSOLVER_MAX_TIME_DER, std::vector<ItemSolverBase*> ())),
//    m_sols ({}),
//    m_dfstore (new DFStore())
//{
//    for (std::size_t i = 0; i < UPSOLVER_MAX_TIME_DER; ++i)
//        m_sols.push_back (new QueueSolutions(&m_store->mesh));
//}

//UpgradableSolverBase::~UpgradableSolverBase ()
//{
//    delete m_dfstore;

//    for (std::vector<ItemSolverBase*> listitem : m_items)
//        for (ItemSolverBase* item : listitem)
//            delete item;

//    for (QueueSolutions* queue : m_sols)
//        delete queue;

//    m_sols.clear ();
//    m_items.clear ();
//}

//void
//UpgradableSolverBase::SetCollectContribution (bool value)
//{
//    m_collectcontributions = value;
//    return;
//}

//void
//UpgradableSolverBase::SetTimeStep (double dt)
//{
//    m_dt = dt;
//    return;
//}

//double
//UpgradableSolverBase::GetTimeStep ()
//{
//    return m_dt;
//}

//void
//UpgradableSolverBase::SetCoeffPen (double coeffPen)
//{
//    m_coeff_pen = coeffPen;
//    return;
//}

//double
//UpgradableSolverBase::GetCoeffPen ()
//{
//    return m_coeff_pen;
//}

//void
//UpgradableSolverBase::SetTime (double time)
//{
//    m_time = time;
//    return;
//}

//double
//UpgradableSolverBase::GetTime ()
//{
//    return m_time;
//}

//void
//UpgradableSolverBase::operator++ ()
//{
//    m_time += m_dt;
//    return;
//}


//std::size_t UpgradableSolverBase::GetMaxActionDerivative ()
//{
//    std::size_t out = 0;

//    for (std::size_t i = 0; i < UPSOLVER_MAX_TIME_DER; ++i)
//        if (this->m_items.at (i).size () != 0)
//            out = i;

//    return out;
//}

//UpgradableSolverBase&
//UpgradableSolverBase::Configure (bool view)
//{
//    INFOS << "Configure" << *this << ENDLINE;

//    for (std::vector<ItemSolverBase*> listitem : m_items)
//        for (ItemSolverBase* item : listitem)
//        {
//            item->Compute (m_time, view, true);


//            if (m_collectcontributions)
//            {
//                std::vector<int>* applicableCells= item->GetApplicableCellsVector ();
//                std::vector<int>* applicableEdges= item->GetApplicableEdgesVector ();

//                std::string name = "Item_"+ToString (item->GetTypeOfItem ())+"_"+to_string (item->GetTagToApply ());

//                m_store->mesh->GetCellsData ()->GetIntArrays ()->Add (name, *applicableCells);

//                m_store->mesh->GetEdgesData ()->GetIntArrays ()->Add (name, *applicableEdges);
//            }
//        }

//    return *this;
//}

//void UpgradableSolverBase::ComputeSecMemberAndMats (std::vector<SparseMatrix>* listMats, PlainVector* secMember, bool view)
//{
//    int numPoints = m_store->mesh->GetNumberOfPoints ();

//    listMats->resize (UPSOLVER_MAX_TIME_DER);
//    for (std::size_t i = 0; i < UPSOLVER_MAX_TIME_DER; ++i)
//        listMats->at (i) = SparseMatrix (numPoints, numPoints);

//    secMember->resize (numPoints);
//    secMember->setZero ();

//    for (std::vector<ItemSolverBase*> listitem : m_items)
//        for (ItemSolverBase* item : listitem)
//        {
//            item->Compute (m_time, view);
//            listMats->at (item->GetTemporalDerivationOrder ())  += *item->CastSparseMatrix ();
//            *secMember += *item->CastPlainVector ();
//        }

//    //    for (std::size_t i = 0; i < listMats->size (); ++i)
//    //    {
//    //        listMats->at (i) = listMats->at (i).pruned (0, 1E-30);
//    //        listMats->at (i).makeCompressed ();
//    //    }

//    return;
//}

//std::ostream& operator<<(std::ostream& out, const UpgradableSolverBase& solver)
//{
//    out << "Solver : "          << COLOR_GREEN  << solver.m_scheme      << FLUSHLINE;
//    out << " [t = "             << COLOR_BLUE   << solver.m_time        << FLUSHLINE;
//    out << ", dt = "            << COLOR_BLUE   << solver.m_dt          << FLUSHLINE;
//    out << ", pen = "           << COLOR_RED    << solver.m_coeff_pen << "/" << solver.m_store->mesh->GetPrescribedSize ()   << COLOR_DEFAULT << "]" << std::flush;

//    return out;
//}
