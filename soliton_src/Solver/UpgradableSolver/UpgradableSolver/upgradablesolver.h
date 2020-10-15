//#ifndef UPGRADABLESOLVERCLASS_H
//#define UPGRADABLESOLVERCLASS_H

//#include <Solver/UpgradableSolver/UpgradableSolverBase/upgradablesolverbase.h>


//template <SCHEME_TEMPORAL_SOLVER t_scheme>
//class UpgradableSolver : public UpgradableSolverBase
//{
//public:
//    UpgradableSolver (Sto4Sol* store);
//    ~UpgradableSolver () override;

//    PlainVector*        Compute() override;
//    UpgradableSolver&   Configure (bool view = true);
//    std::size_t         GetMaxActionDerivative ();

//    template <SCHEME_TEMPORAL_SOLVER t_type>
//    friend std::ostream& operator<<(std::ostream& out, const UpgradableSolver<t_type>& solver);

//private:
//    UpgradableSolver (const UpgradableSolver<t_scheme>& tocopy) = delete;
//    UpgradableSolver& operator= (UpgradableSolver<t_scheme>&) = delete;
//};

//template<SCHEME_TEMPORAL_SOLVER t_scheme>
//UpgradableSolver<t_scheme>::UpgradableSolver (Sto4Sol* store) :
//    UpgradableSolverBase (store)
//{
//    this->m_scheme = t_scheme;
//}

//template<SCHEME_TEMPORAL_SOLVER t_scheme>
//UpgradableSolver<t_scheme>::~UpgradableSolver()
//{}

//UpgradableSolverBase* NewUpgradableSolver (Sto4Sol* store, SCHEME_TEMPORAL_SOLVER scheme);

//template <SCHEME_TEMPORAL_SOLVER t_type>
//std::ostream& operator<<(std::ostream& out, const UpgradableSolver<t_type>& solver)
//{
//    out << "Solver : "          << COLOR_GREEN  << solver.m_scheme      << FLUSHLINE;
//    out << " [t = "             << COLOR_BLUE   << solver.m_time        << FLUSHLINE;
//    out << ", dt = "            << COLOR_BLUE   << solver.m_dt          << FLUSHLINE;
//    out << ", pen = "           << COLOR_RED    << solver.m_coeff_pen << "/" << solver.m_store->mesh->GetPrescribedSize ()   << COLOR_DEFAULT << "]" << std::flush;

//    return out;
//}

//#define EXPLICIT_INST(X)                                           \
//    template<> PlainVector* UpgradableSolver<SCHEME_TEMPORAL_SOLVER::X>::Compute()

//EXPLICIT_INST (NO_TIME);
//EXPLICIT_INST (EULER_EXPLICIT);
//EXPLICIT_INST (EULER_IMPLICIT);
//EXPLICIT_INST (TVD_RK_2);
//EXPLICIT_INST (TVD_RK_4);
//EXPLICIT_INST (NEWMARK_SCHEME);

//#undef EXPLICIT_INST

//#endif // UPGRADABLESOLVERCLASS_H
