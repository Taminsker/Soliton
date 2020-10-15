//#ifndef UPGRADABLESOLVERBASE_H
//#define UPGRADABLESOLVERBASE_H

//#include <tuple>
//#include <vector>
//#include <queue>

//#include <Solver/UpgradableSolver/ItemSolver/itemsolver.h>
//#include <Solver/UpgradableSolver/QueueSolutions/queuesolutions.h>
//#include <Solver/Sto4Sol/sto4sol.h>

//#include <Solver/UpgradableSolver/UpgradableSolverBase/scheme_temporal_solver.h>

//class QueueSolutions;
//class ItemSolverDatStruct;

//class UpgradableSolverBase
//{
//public:
//    UpgradableSolverBase (Sto4Sol* store);
//    virtual ~UpgradableSolverBase ();

//    void    SetCollectContribution (bool value);
//    void    SetTimeStep (double dt);
//    double  GetTimeStep ();
//    void    SetCoeffPen (double coeffPen);
//    double  GetCoeffPen ();
//    void    SetTime (double time);
//    double  GetTime ();
//    void    operator++ ();

//    template<ITEM_SOLVER_TYPE t_type>
//    UpgradableSolverBase&
//    operator<< (ItemSolver<t_type> item);

//    virtual PlainVector*        Compute() = 0;
//    UpgradableSolverBase&       Configure (bool view = true);
//    std::size_t                 GetMaxActionDerivative ();

//    friend std::ostream& operator<<(std::ostream& out, const UpgradableSolverBase& solver);

//protected:
//    /* !!! INTERNAL OBJECTS FOR ITEMS */
//    Sto4Sol*    m_store;
//    double      m_time;
//    double      m_coeff_pen;
//    double      m_dt;
//    bool        m_collectcontributions;

//    /* !!! ITEMS */
//    unstd::matrix<ItemSolverBase*>                m_items;

//    /* !! SOLVERCORE */
//    std::vector<QueueSolutions*>                m_sols;
//    DFStore*                                    m_dfstore;

//    /* !! CONFIGURE */
//    SCHEME_TEMPORAL_SOLVER                      m_scheme;

//    void ComputeSecMemberAndMats (std::vector<SparseMatrix>* listMats,
//                                  PlainVector* secMember, bool view = false);

//private:
//    UpgradableSolverBase (const UpgradableSolverBase& tocopy) = delete;
//    UpgradableSolverBase& operator= (UpgradableSolverBase&) = delete;
//};

//template<ITEM_SOLVER_TYPE t_type>
//UpgradableSolverBase&
//UpgradableSolverBase::operator<< (ItemSolver<t_type> item)
//{
//    ItemSolver<t_type>* newitem = new ItemSolver<t_type>(item);

//    newitem->ConfigureFromSolver (&m_store->mesh, &m_coeff_pen, &m_collectcontributions);

//    std::size_t actionDerivative = newitem->GetTemporalDerivationOrder ();
//    m_items.at (actionDerivative).push_back (newitem);

//    return *this;
//}

//std::ostream& operator<<(std::ostream& out, const UpgradableSolverBase& solver);


//#endif // UPGRADABLESOLVERBASE_H
