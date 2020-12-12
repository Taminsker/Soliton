#ifndef SRC_SOLVER_SOLITONFE_SOLITONFEITEMBASE_SOLITONFEITEMBASE_HPP
#define SRC_SOLVER_SOLITONFE_SOLITONFEITEMBASE_SOLITONFEITEMBASE_HPP

#include "../../../Algorithms/Math/math.hpp"
#include "../../../Core/core.hpp"
#include "../../../Enums/enums.hpp"
#include "../../DFStruct/dfstruct.hpp"
#include "../../FEStruct/festruct.hpp"
#include "../../MeshStorage/meshstorage.hpp"
#include "../../QuadStruct/quadstruct.hpp"
#include "../../SolitonFunctor/solitonfunctor.hpp"

#define NULL_FUNC_NAME soliton_null_fun
real_t NULL_FUNC_NAME (Point a, real_t b);

class SolitonFEContainer;

class SolitonFEItemBase
{
public:
    explicit SolitonFEItemBase (SolitonFEContainer * parent, PHYS tag = PHYS::DEFAULT);
    SolitonFEItemBase (const SolitonFEItemBase & tocopy);
    virtual ~SolitonFEItemBase ();

    virtual void   SetAdditionalName (std::string name) = 0;
    virtual ITEM_T GetType () const                     = 0;
    virtual void   Compute (real_t time)                = 0;

    void SetSparseMatrix (std::vector<Triplet> * listTriplets);
    void SetSecMember (std::vector<Duo> * listDuos);

    SOLITON_INLINE
    void
    SetFunction2Apply (SolitonFunctor * fun)
    {
        m_fun = fun;
        return;
    }

    SOLITON_INLINE
    void
    SetVariableOverTimeOn ()
    {
        m_var0vrT = true;
        return;
    }

    SOLITON_INLINE
    void
    SetVariableOverTimeOff ()
    {
        m_var0vrT = true;
        return;
    }

    SOLITON_INLINE
    void
    SetVariableOverTime (bool var)
    {
        m_var0vrT = var;
        return;
    }

    SOLITON_INLINE
    void
    SetTemporalOrderDerivative (ul_t t)
    {
        m_derT = t;
        return;
    }

    SOLITON_INLINE
    ul_t
    GetTemporalOrderDerivative () const
    {
        return m_derT;
    }

    SOLITON_INLINE
    void
    SetIdTarget (ul_t id)
    {
        m_id_target = id;
        return;
    }

    SOLITON_INLINE
    const SparseMatrix *
    GetSparseMatrix () const
    {
        return m_mat;
    }

    SOLITON_INLINE
    const DenseVector *
    GetPlainVector () const
    {
        return m_sec;
    }

    SOLITON_INLINE
    std::string
    GetName () const
    {
        return m_name;
    }

    SOLITON_INLINE
    SolitonFEContainer *
    GetContainer () const
    {
        return m_c;
    }

    void
    Print (std::ostream & out);

protected:
    SolitonFEContainer * m_c;
    ul_t                 m_id_target;

    SparseMatrix *    m_mat;
    DenseVector *     m_sec;
    std::vector<int> *m_appCells, *m_appEdges, *m_appPoints;

    std::string      m_name;
    PHYS             m_tag2Apply;
    bool             m_var0vrT;
    ul_t             m_derT;
    SolitonFunctor * m_fun;
};

#endif /* SRC_SOLVER_SOLITONFE_SOLITONFEITEMBASE_SOLITONFEITEMBASE_HPP */
