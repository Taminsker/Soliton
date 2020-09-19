#ifndef UPGRADABLESOLVERCLASS_H
#define UPGRADABLESOLVERCLASS_H

#include <tuple>
#include <vector>
#include <queue>

#include <Solver/UpgradableSolver/ItemSolver/itemsolver.h>
#include <Solver/UpgradableSolver/QueueSolutions/queuesolutions.h>
#include <Solver/Sto4Sol/sto4sol.h>

#include "scheme_temporal_solver.h"


class QueueSolutions;

template <SCHEME_TEMPORAL_SOLVER t_scheme>
class UpgradableSolver
{
public:
    UpgradableSolver (Sto4Sol* store);
    ~UpgradableSolver ();

    void    SetTimeStep (double dt);
    double  GetTimeStep ();
    void    SetCoeffPen (double coeffPen);
    double  GetCoeffPen ();
    void    SetTime (double time);
    double  GetTime ();
    void    operator++ ();

    template<ITEM_SOLVER_TYPE t_type>
    UpgradableSolver&
    operator<< (ItemSolver<t_type> item);

    PlainVector*        Compute();
    UpgradableSolver&   Configure (bool view);
    std::size_t         GetMaxActionDerivative ();


    template <SCHEME_TEMPORAL_SOLVER t_type>
    friend std::ostream& operator<<(std::ostream& out, const UpgradableSolver<t_type>& solver);

protected:
    /* !!! INTERNAL OBJECTS FOR ITEMS */
    Sto4Sol*    m_store;
    Cell        *m_current_cell;
    double      m_time;
    double      m_coeff_pen;
    double      m_dt;

    /* !!! ITEMS */
    std::matrix<ItemSolverBase*>                m_items;

    /* !! SOLVERCORE */
    std::vector<QueueSolutions*>                m_sols;
    DFStore*                                    m_dfstore;

    /* !! CONFIGURE */
    SCHEME_TEMPORAL_SOLVER                      m_scheme;

    void ComputeSecMemberAndMats (std::vector<SparseMatrix>* listMats, PlainVector* secMember, bool view = false);

private:
    UpgradableSolver (const UpgradableSolver<t_scheme>& tocopy) = delete;
    UpgradableSolver& operator= (UpgradableSolver<t_scheme>&) = delete;
};

template<SCHEME_TEMPORAL_SOLVER t_scheme>
UpgradableSolver<t_scheme>::UpgradableSolver (Sto4Sol* store) :
    m_store (store),
    m_current_cell (nullptr),
    m_time (0.),
    m_coeff_pen (1.),
    m_dt (1.),
    m_items (std::matrix<ItemSolverBase*> (UPSOLVER_MAX_TIME_DER, std::vector<ItemSolverBase*> ())),
    m_sols ({}),
    m_dfstore (new DFStore()),
    m_scheme (t_scheme)
{
    for (std::size_t i = 0; i < UPSOLVER_MAX_TIME_DER; ++i)
        m_sols.push_back (new QueueSolutions(&m_store->mesh));
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
UpgradableSolver<t_scheme>::~UpgradableSolver()
{
    delete m_dfstore;

    for (std::vector<ItemSolverBase*> listitem : m_items)
        for (ItemSolverBase* item : listitem)
            delete item;

    for (QueueSolutions* queue : m_sols)
        delete queue;

    m_sols.clear ();
    m_items.clear ();
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
void
UpgradableSolver<t_scheme>::SetTimeStep (double dt)
{
    m_dt = dt;
    return;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
double
UpgradableSolver<t_scheme>::GetTimeStep ()
{
    return m_dt;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
void
UpgradableSolver<t_scheme>::SetCoeffPen (double coeffPen)
{
    m_coeff_pen = coeffPen;
    return;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
double
UpgradableSolver<t_scheme>::GetCoeffPen ()
{
    return m_coeff_pen;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
void
UpgradableSolver<t_scheme>::SetTime (double time)
{
    m_time = time;
    return;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
double
UpgradableSolver<t_scheme>::GetTime ()
{
    return m_time;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
void
UpgradableSolver<t_scheme>::operator++ ()
{
    m_time += m_dt;
    return;
}


template<SCHEME_TEMPORAL_SOLVER t_scheme>
std::size_t UpgradableSolver<t_scheme>::GetMaxActionDerivative ()
{
    std::size_t out = 0;

    for (std::size_t i = 0; i < UPSOLVER_MAX_TIME_DER; ++i)
        if (this->m_items.at (i).size () != 0)
            out = i;

    return out;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
UpgradableSolver<t_scheme>&
UpgradableSolver<t_scheme>::Configure (bool view)
{
    INFOS << "Configure" << *this << ENDLINE;

    for (std::vector<ItemSolverBase*> listitem : m_items)
        for (ItemSolverBase* item : listitem)
            item->Compute (m_time, view, true);

    return *this;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
template<ITEM_SOLVER_TYPE t_type>
UpgradableSolver<t_scheme>&
UpgradableSolver<t_scheme>::operator<< (ItemSolver<t_type> item)
{
    ItemSolver<t_type>* newitem = new ItemSolver<t_type>(item);

    newitem->ConfigureFromSolver (&m_store->mesh, &m_coeff_pen);

    std::size_t actionDerivative = newitem->GetTemporalDerivationOrder ();
    m_items.at (actionDerivative).push_back (newitem);

    return *this;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
void UpgradableSolver<t_scheme>::ComputeSecMemberAndMats (std::vector<SparseMatrix>* listMats, PlainVector* secMember, bool view)
{
    int numPoints = m_store->mesh->GetNumberOfPoints ();

    listMats->resize (UPSOLVER_MAX_TIME_DER);
    for (std::size_t i = 0; i < UPSOLVER_MAX_TIME_DER; ++i)
        listMats->at (i) = SparseMatrix (numPoints, numPoints);

    secMember->resize (numPoints);
    secMember->setZero ();

    for (std::vector<ItemSolverBase*> listitem : m_items)
        for (ItemSolverBase* item : listitem)
        {
            item->Compute (m_time, view);

            listMats->at (item->GetTemporalDerivationOrder ())  += *item->CastSparseMatrix ();
            *secMember += *item->CastPlainVector ();
        }

    return;
}

template<SCHEME_TEMPORAL_SOLVER t_scheme>
std::ostream& operator<<(std::ostream& out, const UpgradableSolver<t_scheme>& solver)
{
    out << "Solver : "          << COLOR_GREEN  << solver.m_scheme      << FLUSHLINE;
    out << " [t = "             << COLOR_BLUE   << solver.m_time        << FLUSHLINE;
    out << ", dt = "            << COLOR_BLUE   << solver.m_dt          << FLUSHLINE;
    out << ", pen = "           << COLOR_RED    << solver.m_coeff_pen   << COLOR_DEFAULT << "]" << std::flush;

    return out;
}


#define EXPLICIT_INST(X)                                           \
    template<> PlainVector* UpgradableSolver<SCHEME_TEMPORAL_SOLVER::X>::Compute()

EXPLICIT_INST (NO_TIME);
EXPLICIT_INST (EULER_EXPLICIT);
EXPLICIT_INST (EULER_IMPLICIT);
EXPLICIT_INST (TVD_RK_2);
EXPLICIT_INST (TVD_RK_4);
EXPLICIT_INST (NEWMARK_SCHEME);

#undef EXPLICIT_INST

#endif // UPGRADABLESOLVERCLASS_H
