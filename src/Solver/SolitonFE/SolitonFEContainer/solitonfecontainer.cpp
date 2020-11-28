#include "solitonfecontainer.h"

SolitonFEContainer::SolitonFEContainer (MeshStorage *store, INTER tag) :
    m_tag_intersection (tag),
    m_store (store),
    m_time (0x0),
    m_coeff_pen (1E2),
    m_pow_pen (1),
    m_dt (0x1),
    m_collect_contributions (false),
    m_list_items (new unstd::matrix<SolitonFEItemBase*> (SOLITONCONTAINER_ORDER_MAX)),
    m_queues ({}),
    m_dfstore (new DFStore ())
{
    for (ul_t i = 0; i < SOLITONCONTAINER_ORDER_MAX; ++i)
        m_queues.push_back (new SolitonQueue(m_store->GetMainMesh_inc ()));
}

SolitonFEContainer::~SolitonFEContainer ()
{
    for (auto li : *m_list_items)
        for (auto item : li)
            delete item;

    delete m_list_items;

    for (auto q : m_queues)
        delete q;

    m_queues.clear ();

    delete m_dfstore;
}
