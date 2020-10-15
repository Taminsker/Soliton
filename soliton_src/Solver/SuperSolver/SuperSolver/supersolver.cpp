#include "supersolver.h"

SuperSolver::SuperSolver (Sto4Sol *store, std::vector<INTER> types) :
    m_nextsolver (nullptr),
    m_type (SCH_T::DEFAULT),
    m_store (store),
    m_time (0x0),
    m_coeff_pen (1E2),
    m_dt (0x1),
    m_collect_contributions (false),
    m_list_items (new unstd::matrix<SuperItem*> (SUPERSOLVER_MAX)),
    m_queue (new SuperQueue(&m_store->mesh)),
    m_dfstore (new DFStore ())
{
    if (types.empty ())
    {
        INFOS << "Default tag intersection in solver is : " << to_string(INTER::DEFAULT) << ENDLINE;
        return;
    }

    m_tag_intersection = types.front ();

    if (types.size () > 1)
        m_nextsolver = new SuperSolver (store, std::vector<INTER>(types.begin () + 0x1, types.end ()));
}

SuperSolver::SuperSolver (const SuperSolver& tocopy) :
    m_nextsolver (nullptr),
    m_tag_intersection (tocopy.m_tag_intersection),
    m_type (tocopy.m_type),
    m_store (tocopy.m_store),
    m_time (tocopy.m_time),
    m_coeff_pen (tocopy.m_coeff_pen),
    m_dt (tocopy.m_dt),
    m_collect_contributions (tocopy.m_collect_contributions),
    m_list_items (new unstd::matrix<SuperItem*> (SUPERSOLVER_MAX)),
    m_queue (new SuperQueue(&m_store->mesh)),
    m_dfstore (new DFStore ())
{
    if (tocopy.m_nextsolver != nullptr)
        m_nextsolver = new SuperSolver (*tocopy.m_nextsolver);

    for (std::vector<SuperItem*> listitem : *tocopy.m_list_items)
        for (SuperItem* item : listitem)
            AddItem (*item, item->GetOrder ());
}

SuperSolver::~SuperSolver ()
{
    delete m_list_items;
    delete m_queue;
    delete m_dfstore;

    if (m_nextsolver != nullptr)
        delete m_nextsolver;
}
