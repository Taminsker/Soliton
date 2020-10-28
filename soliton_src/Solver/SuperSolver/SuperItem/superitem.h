#ifndef SUPERITEM_H
#define SUPERITEM_H

#include <Core/core.h>
#include <Enums/enums.h>
#include <Solver/DFStruct/dfstruct.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>

#define NULL_FUNC_NAME soliton_null_fun
double NULL_FUNC_NAME (Point a, double b);

class SuperSolver;
class SuperItem;

//template<ITEM_T t_type>
//void Compute (SuperItem* item);

class SuperItem
{
public:   
    SuperItem (ITEM_T type);
    SuperItem (const SuperItem& tocopy);
    ~SuperItem ();

    void SetSuperSolver (SuperSolver *solver);

    SOLITON_INLINE
    void SetTag2Apply (PHYS tag)
    {
        m_tag2Apply = tag;
        return;
    }

    SOLITON_INLINE
    void SetFunction2Apply (std::function<double(Point, double)> fun)
    {
        m_fun = fun;
        return;
    }

    SOLITON_INLINE
    void SetDerT (std::size_t derT)
    {
        m_derT = derT;
        return;
    }

    SOLITON_INLINE
    void SetVariableOverTime (bool variableOverTime)
    {
        m_var0vrT = variableOverTime;
        return;
    }

    SOLITON_INLINE
    std::size_t GetOrder () const
    {
        return m_derT;
    }

    SOLITON_INLINE
    bool ItsVariableOverTime () const
    {
        return m_var0vrT;
    }

    SOLITON_INLINE
    PHYS GetTag2Apply() const
    {
        return m_tag2Apply;
    }

    SOLITON_INLINE
    const SparseMatrix* GetSparseMatrix () const
    {
        return m_mat;
    }

    SOLITON_INLINE
    const PlainVector* GetPlainVector () const
    {
        return m_sec;
    }

    SOLITON_INLINE
    void SetSparseMatrix (std::vector<Triplet>* listTriplets)
    {
        int numPoints  = (*m_mesh)->GetNumberOfPoints ();
        m_mat->resize (numPoints, numPoints);
        m_mat->setFromTriplets(listTriplets->begin (), listTriplets->end ());
        return;
    }

    SOLITON_INLINE
    void SetSecMember (std::vector<Duo>* listDuos)
    {
        int numPoints  = (*m_mesh)->GetNumberOfPoints ();

        m_sec->setZero (numPoints);

        for (Duo duo : *listDuos)
            m_sec->coeffRef (duo.idx ()) += duo.value ();

        return;
    }

    void Compute ();

    SOLITON_INLINE
    ITEM_T GetTypeOfItem() const
    {
        return m_type;
    }

    SOLITON_INLINE
    Mesh* GetMesh () const
    {
        return *m_mesh;
    }

    SOLITON_INLINE
    std::function<double(Point, double)>
    GetFun ()
    {
        return m_fun;
    }

    SOLITON_INLINE
    double
    GetTime ()
    {
        return *m_time;
    }

    SOLITON_INLINE
    bool
    GetView ()
    {
        return *m_view;
    }

    SOLITON_INLINE
    bool
    GetCollect ()
    {
        return *m_collect;
    }

    SOLITON_INLINE
    bool
    GetForce ()
    {
        return *m_force_compute;
    }

    SOLITON_INLINE
    double
    GetCoeffPen ()
    {
        return *m_coeffpen;
    }

    SOLITON_INLINE
    int GetPowPen ()
    {
        return *m_powpen;
    }

    SOLITON_INLINE
    void
    SetTargetMesh (Mesh **target)
    {
        m_target = target;
        return;
    }

    SOLITON_INLINE
    Mesh*
    GetTargetMesh ()
    {
        if (m_target != nullptr)
            return *m_target;
        return nullptr;
    }

    SOLITON_INLINE
    SuperSolver*
    GetSolver ()
    {
        return m_solver;
    }

    SOLITON_INLINE
    std::vector<int>* GetApplicationCell ()
    {
        return m_appCells;
    }

    SOLITON_INLINE
    std::vector<int>* GetApplicationEdge ()
    {
        return m_appEdges;
    }

    SOLITON_INLINE
    std::vector<int>* GetApplicationPoint ()
    {
        return m_appPoints;
    }

    SOLITON_INLINE
    std::string GetWholeName () const
    {
        return m_name;
    }

    SOLITON_INLINE
    void SetAdditionalName (std::string name)
    {
        m_name = to_string(m_type) + name;
        return;
    }

    friend std::ostream& operator<< (std::ostream& out, const SuperItem& item);

    friend class SuperItemCompute;

private:

    // INTERNAL COMPUTE
    template<ITEM_T t_type>
    void InternalCompute ()
    {}

protected:
    const ITEM_T                            m_type;
    Mesh                                    **m_mesh;
    Mesh                                    **m_target;
    SuperSolver                             *m_solver;
    SparseMatrix                            *m_mat;
    PlainVector                             *m_sec;
    std::vector<int>                        *m_appCells, *m_appEdges, *m_appPoints;

//    std::size_t                             m_id_item;
    std::string                             m_name;
    PHYS                                    m_tag2Apply;
    bool                                    m_var0vrT, *m_collect, *m_view, *m_force_compute;
    double                                  *m_coeffpen, *m_time;
    int                                     *m_powpen;
    std::size_t                             m_derT;
    std::function<double(Point, double)>    m_fun;
};

#define EXPLICIT_INST(X) \
    template<> void SuperItem::InternalCompute<ITEM_T::X> ()

EXPLICIT_INST(EMPTY);
EXPLICIT_INST(PHI_PHI);
EXPLICIT_INST(DAMPING_PHI_PHI);
EXPLICIT_INST(GRAD_GRAD);
EXPLICIT_INST(SECOND_MEMBER);
EXPLICIT_INST(DIRICHLET_BOUNDARY);
EXPLICIT_INST(NEUMANN_BOUNDARY);
EXPLICIT_INST(DIRICHLET_BOUNDARY_SBM);
EXPLICIT_INST(NEUMANN_BOUNDARY_SBM);

#undef EXPLICIT_INST

SOLITON_INLINE
void SuperItem::Compute()
{
    switch (m_type)
    {
    case ITEM_T::EMPTY:
        return this->InternalCompute<ITEM_T::EMPTY> ();
    case ITEM_T::PHI_PHI:
        return this->InternalCompute<ITEM_T::PHI_PHI> ();
    case ITEM_T::DAMPING_PHI_PHI:
        return this->InternalCompute<ITEM_T::DAMPING_PHI_PHI> ();
    case ITEM_T::GRAD_GRAD:
        return this->InternalCompute<ITEM_T::GRAD_GRAD> ();
    case ITEM_T::SECOND_MEMBER:
        return this->InternalCompute<ITEM_T::SECOND_MEMBER> ();
    case ITEM_T::DIRICHLET_BOUNDARY:
        return this->InternalCompute<ITEM_T::DIRICHLET_BOUNDARY> ();
    case ITEM_T::NEUMANN_BOUNDARY:
        return this->InternalCompute<ITEM_T::NEUMANN_BOUNDARY> ();
    case ITEM_T::DIRICHLET_BOUNDARY_SBM:
        return this->InternalCompute<ITEM_T::DIRICHLET_BOUNDARY_SBM> ();
    case ITEM_T::NEUMANN_BOUNDARY_SBM:
        return this->InternalCompute<ITEM_T::NEUMANN_BOUNDARY_SBM> ();
    }

    ERROR << "Type incorrect ??" << ENDLINE;
    return;
}
#endif // SUPERITEM_H
