#ifndef SRC_SOLVER_SOLITONFE_SOLITONFECONTAINER_SOLITONFECONTAINER_HPP
#define SRC_SOLVER_SOLITONFE_SOLITONFECONTAINER_SOLITONFECONTAINER_HPP

#include "../../../Algorithms/Math/math.hpp"
#include "../../../Enums/enums.hpp"
#include "../../DFStruct/dfstruct.hpp"
#include "../../MeshStorage/meshstorage.hpp"
#include "../../SolitonQueue/solitonqueue.hpp"
#include "../SolitonFEItem/solitonfeitem.hpp"

class SolitonFEContainer
{
public:
    SolitonFEContainer (MeshStorage * store, INTER tag);
    ~SolitonFEContainer ();

    SOLITON_INLINE
    INTER
    GetInterTag ()
    {
        return m_tag_intersection;
    }

    SOLITON_INLINE
    void
    SetCollectContribution (bool value)
    {
        m_collect_contributions = value;
        return;
    }

    SOLITON_INLINE
    void
    SetCollectContributionOn ()
    {
        m_collect_contributions = true;
        return;
    }

    SOLITON_INLINE
    void
    SetCollectContributionOff ()
    {
        m_collect_contributions = false;
        return;
    }

    SOLITON_INLINE
    void
    SetForceCompute (bool force)
    {
        m_force_compute = force;
        return;
    }

    SOLITON_INLINE
    void
    SetForceComputeOn ()
    {
        m_force_compute = true;
        return;
    }

    SOLITON_INLINE
    void
    SetForceComputeOff ()
    {
        m_force_compute = false;
        return;
    }

    SOLITON_INLINE
    void
    SetView (bool view)
    {
        m_view = view;
        return;
    }

    SOLITON_INLINE
    void
    SetViewOn ()
    {
        m_view = true;
        return;
    }

    SOLITON_INLINE
    void
    SetViewOff ()
    {
        m_view = false;
        return;
    }

    SOLITON_INLINE
    void
    SetCoeffPen (real_t coeffPen)
    {
        m_coeff_pen = coeffPen;
        return;
    }

    SOLITON_INLINE
    real_t
    GetCoeffPen ()
    {
        return m_coeff_pen;
    }

    SOLITON_INLINE
    void
    SetPowPen (int powPen)
    {
        m_pow_pen = powPen;
        return;
    }

    SOLITON_INLINE
    real_t
    GetPowPen ()
    {
        return m_pow_pen;
    }

    SOLITON_INLINE
    void
    Print ()
    {
        INFOS << "Solver : " << COLOR_BLUE << to_string (m_tag_intersection) << COLOR_DEFAULT << FLUSHLINE;
        COUT << " [t = " << COLOR_RED << m_time << FLUSHLINE;
        COUT << ", dt = " << COLOR_BLUE << m_dt << FLUSHLINE;
        COUT << ", pen = " << COLOR_RED << m_coeff_pen << "/" << m_store->GetMainMesh ()->GetPrescribedSize () << "^" << m_pow_pen << COLOR_DEFAULT << FLUSHLINE;
        COUT << ", numItems = " << COLOR_BLUE << GetNumberOfItems () << COLOR_DEFAULT << "]" << FLUSHLINE;
        COUT << ENDLINE;

        return;
    }

    SOLITON_INLINE
    void
    ComputeSecMemberAndMats (std::vector<SparseMatrix> * listMats, DenseVector * secMember, real_t time)
    {
        int numPoints = m_store->GetMainMesh ()->GetNumberOfPoints ();

        if (m_view)
            this->Print ();

        listMats->resize (SOLITONCONTAINER_ORDER_MAX);
        for (uid_t i = 0; i < SOLITONCONTAINER_ORDER_MAX; ++i)
            listMats->at (i) = SparseMatrix (numPoints, numPoints);

        secMember->resize (numPoints);
        secMember->setZero ();

        for (std::vector<SolitonFEItemBase *> listitem : *m_list_items)
            for (SolitonFEItemBase * item : listitem)
            {
                item->Compute (time);
                listMats->at (item->GetTemporalOrderDerivative ()) += *item->GetSparseMatrix ();
                *secMember += *item->GetPlainVector ();
            }

        return;
    }

    template <ITEM_T t_type>
    SOLITON_INLINE void
    AddItem (SolitonFEItem<t_type> item)
    {
        ul_t alpha = item.GetTemporalOrderDerivative ();

        SolitonFEItem<t_type> * toadd = new SolitonFEItem<t_type> (item);
        toadd->SetAdditionalName ("_" + std::to_string (alpha * SOLITONCONTAINER_ORDER_MAX + m_list_items->at (alpha).size ()) + "_" + to_string (m_tag_intersection));
        m_list_items->at (alpha).push_back (toadd);

        return;
    }

    SOLITON_INLINE
    int
    GetNumberOfItems ()
    {
        int count = 0;

        for (auto lt : *m_list_items)
            count += lt.size ();

        return count;
    }

    SOLITON_INLINE
    void
    Configure (real_t time = 0x0)
    {
        if (m_view)
        {
            BEGIN << "Configure solver" << ENDLINE;
            this->Print ();
        }

        for (std::vector<SolitonFEItemBase *> listitem : *m_list_items)
            for (SolitonFEItemBase * item : listitem)
                item->Compute (time);

        ENDFUN;
        return;
    }

    SOLITON_INLINE
    SolitonQueue *
    GetQueue (ul_t id)
    {
        return m_queues.at (id);
    }

    friend class SolitonFEItemBase;

    template <ITEM_T U>
    friend class SolitonFEItem;

protected:
    INTER         m_tag_intersection;
    MeshStorage * m_store;
    real_t        m_time;
    real_t        m_coeff_pen;
    int           m_pow_pen;
    real_t        m_dt;
    bool          m_collect_contributions, m_view, m_force_compute;

    unstd::matrix<SolitonFEItemBase *> * m_list_items;
    std::vector<SolitonQueue *>          m_queues;
    DFStore *                            m_dfstore;

private:
    SolitonFEContainer (const SolitonFEContainer &) = delete;
};

#endif /* SRC_SOLVER_SOLITONFE_SOLITONFECONTAINER_SOLITONFECONTAINER_HPP */
