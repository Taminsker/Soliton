#ifndef SUPERSOLVERCLASS_H
#define SUPERSOLVERCLASS_H

#include <iostream>

#include <Algorithms/Math/math.h>
#include <Solver/DFStruct/dfstruct.h>
#include <Solver/Sto4Sol/sto4sol.h>

#include <Solver/SuperSolver/SuperItem/superitem.h>
#include <Solver/SuperSolver/SuperQueue/superqueue.h>

class SuperSolver
{
public:
    SuperSolver (Sto4Sol *store, std::vector<INTER> types);
    SuperSolver (const SuperSolver& tocopy);
    ~SuperSolver ();

    SOLITON_INLINE
    INTER
    GetInterTag ()
    {
        return m_tag_intersection;
    }

    SOLITON_INLINE
    std::vector<INTER>
    GetListOfInterTags ()
    {
        std::vector<INTER> out = {};
        if (m_nextsolver != nullptr)
            out = m_nextsolver->GetListOfInterTags ();

        out.push_back (GetInterTag ());
        return out;
    }

    SOLITON_INLINE
    int GetDepth ()
    {
        if (m_nextsolver != nullptr)
            return 1 + m_nextsolver->GetDepth ();
        return 1;
    }

    SOLITON_INLINE
    void SetCollectContribution (bool value)
    {
        m_collect_contributions = value;
        if (m_nextsolver != nullptr)
            return m_nextsolver->SetCollectContribution (value);
        return;
    }

    SOLITON_INLINE
    void SetCollectContributionOn ()
    {
        m_collect_contributions = true;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetCollectContributionOn ();
        return;
    }

    SOLITON_INLINE
    void SetCollectContributionOff ()
    {
        m_collect_contributions = false;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetCollectContributionOff ();
        return;
    }

    SOLITON_INLINE
    void SetForceCompute (bool force)
    {
        m_force_compute = force;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetForceCompute (force);
        return;
    }

    SOLITON_INLINE
    void SetForceComputeOn ()
    {
        m_force_compute = true;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetForceComputeOn ();
        return;
    }

    SOLITON_INLINE
    void SetForceComputeOff ()
    {
        m_force_compute = false;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetForceComputeOff ();
        return;
    }

    SOLITON_INLINE
    void SetView (bool view)
    {
        m_view = view;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetView (view);
        return;
    }

    SOLITON_INLINE
    void SetViewOn ()
    {
        m_view = true;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetViewOn ();
        return;
    }

    SOLITON_INLINE
    void SetViewOff ()
    {
        m_view = false;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetViewOff ();
        return;
    }

    SOLITON_INLINE
    void SetTimeStep (double dt)
    {
        m_dt = dt;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetTimeStep (dt);
        return;
    }

    SOLITON_INLINE
    double GetTimeStep ()
    {
        return m_dt;
    }

    SOLITON_INLINE
    void SetCoeffPen (double coeffPen)
    {
        m_coeff_pen = coeffPen;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetCoeffPen (coeffPen);
        return;
    }

    SOLITON_INLINE
    double GetCoeffPen ()
    {
        return m_coeff_pen;
    }

    SOLITON_INLINE
    void SetPowPen (int powPen)
    {
        m_pow_pen = powPen;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetPowPen (powPen);
        return;
    }

    SOLITON_INLINE
    double GetPowPen ()
    {
        return m_pow_pen;
    }

    SOLITON_INLINE
    void SetTime (double time)
    {
        m_time = time;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetTime (time);
        return;
    }

    SOLITON_INLINE
    double GetTime ()
    {
        return m_time;
    }

    SOLITON_INLINE
    void operator++ ()
    {
        m_time += m_dt;

        if (m_nextsolver != nullptr)
            return m_nextsolver->operator++ ();
        return;
    }

    SOLITON_INLINE
    void operator+= (double dt)
    {
        m_time += dt;

        if (m_nextsolver != nullptr)
            return m_nextsolver->operator+= (dt);
        return;
    }

    SOLITON_INLINE
    void Print ()
    {
        INFOS << "Solver : id-depth " << COLOR_RED << GetDepth () << COLOR_DEFAULT << " : " << FLUSHLINE;
        COUT << COLOR_BLUE << to_string(m_tag_intersection) << COLOR_DEFAULT << FLUSHLINE;
        COUT << " [t = " << COLOR_RED << m_time << FLUSHLINE;
        COUT << ", dt = " << COLOR_BLUE << m_dt << FLUSHLINE;
        COUT << ", pen = " << COLOR_RED << m_coeff_pen << "/" << m_store->mesh->GetPrescribedSize () << "^" << m_pow_pen << COLOR_DEFAULT << FLUSHLINE;
        COUT << ", numItems = " << COLOR_BLUE << GetNumberOfItems () << COLOR_DEFAULT <<"]" << FLUSHLINE;
        COUT << ENDLINE;

        return;
    }

    SOLITON_INLINE
    void ComputeSecMemberAndMats (std::vector<SparseMatrix>* listMats, PlainVector* secMember)
    {
        int numPoints = m_store->mesh->GetNumberOfPoints ();

        if (m_view)
            this->Print ();

        listMats->resize (SUPERSOLVER_MAX);
        for (std::size_t i = 0; i < SUPERSOLVER_MAX; ++i)
            listMats->at (i) = SparseMatrix (numPoints, numPoints);

        secMember->resize (numPoints);
        secMember->setZero ();

        for (std::vector<SuperItem*> listitem : *m_list_items)
            for (SuperItem* item : listitem)
            {
                item->Compute ();
                listMats->at (item->GetOrder ()) += *item->GetSparseMatrix ();
                *secMember += *item->GetPlainVector ();
            }

        return;
    }

    SOLITON_INLINE
    void AddItem (SuperItem item, std::size_t alpha, INTER tag)
    {
        if (tag == m_tag_intersection || tag == INTER::DEFAULT)
        {
            SuperItem* toadd = new SuperItem (item);
            toadd->SetSuperSolver (this);
            toadd->SetAdditionalName ("_" + std::to_string (alpha * SUPERSOLVER_MAX + m_list_items->at (alpha).size ()) + "_" + to_string(m_tag_intersection));
            m_list_items->at (alpha).push_back (toadd);
        }

        if (m_nextsolver != nullptr)
            return m_nextsolver->AddItem (item, alpha, tag);
        return;
    }

    SOLITON_INLINE
    void SetScheme (SCH_T type)
    {
        m_type = type;

        if (m_nextsolver != nullptr)
            return m_nextsolver->SetScheme (type);
        return;
    }

    SOLITON_INLINE
    SCH_T GetScheme ()
    {
        return m_type;
    }

    SOLITON_INLINE
    int GetNumberOfItems ()
    {
        int count = 0;

        for (auto lt : *m_list_items)
            count += lt.size ();

        return count;
    }

    SOLITON_INLINE
    void Configure ()
    {
        BEGIN << "Configure solver" << ENDLINE;
        if (m_view)
            this->Print ();

        for (std::vector<SuperItem*> listitem : *m_list_items)
            for (SuperItem* item : listitem)
                item->Compute ();

        ENDFUN;

        if (m_nextsolver != nullptr)
            return m_nextsolver->Configure ();
        return;
    }

    std::vector<PlainVector*>
    Compute ();

    friend class SuperItem;

private:
    template<SCH_T>
    void InternalCompute ();

protected:
    SuperSolver                     *m_nextsolver;

    INTER                           m_tag_intersection;
    SCH_T                           m_type;
    Sto4Sol*                        m_store;
    double                          m_time;
    double                          m_coeff_pen;
    int                             m_pow_pen;
    double                          m_dt;
    bool                            m_collect_contributions, m_view, m_force_compute;

    unstd::matrix<SuperItem*>       *m_list_items;
    std::vector<SuperQueue*>         m_queues;
    DFStore                         *m_dfstore;
};

#define EXPLICIT_INST(X) \
    template<> void SuperSolver::InternalCompute<SCH_T::X> ()

EXPLICIT_INST(NO_TIME);
EXPLICIT_INST(TVD_RK_2);
EXPLICIT_INST(TVD_RK_4);
EXPLICIT_INST(EULER_EXPLICIT);
EXPLICIT_INST(EULER_IMPLICIT);
EXPLICIT_INST(NEWMARK_SCHEME);

#undef EXPLICIT_INST

SOLITON_INLINE
std::vector<PlainVector*>
SuperSolver::Compute ()
{
    BEGIN << "Compute solver with scheme " << COLOR_GREEN << to_string(this->m_type) << " on " << to_string(this->m_tag_intersection)<< ENDLINE;

    switch (m_type)
    {
    case SCH_T::NO_TIME:
        InternalCompute<SCH_T::NO_TIME> ();
        break;
    case SCH_T::TVD_RK_2:
        InternalCompute<SCH_T::TVD_RK_2> ();
        break;
    case SCH_T::TVD_RK_4:
        InternalCompute<SCH_T::TVD_RK_4> ();
        break;
    case SCH_T::EULER_EXPLICIT:
        InternalCompute<SCH_T::EULER_EXPLICIT> ();
        break;
    case SCH_T::EULER_IMPLICIT:
        InternalCompute<SCH_T::EULER_IMPLICIT> ();
        break;
    case SCH_T::NEWMARK_SCHEME:
        InternalCompute<SCH_T::NEWMARK_SCHEME> ();
        break;
    }

    ENDFUN;

    std::vector<PlainVector*> out = {};
    if (m_nextsolver != nullptr)
        out = m_nextsolver->Compute ();

    out.push_back (m_queues.at (0)->GetNumber (0));
    return out;
}

#endif // SUPERSOLVERCLASS_H
