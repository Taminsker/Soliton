//#ifndef ITEMSOLVER_H
//#define ITEMSOLVER_H

//#include <Solver/UpgradableSolver/ItemSolverBase/itemsolverbase.h>

//template <ITEM_SOLVER_TYPE t_type>
//class ItemSolver : public ItemSolverBase
//{
//public:
//    ItemSolver ();
//    ItemSolver (const ItemSolver& tocopy);
//    ~ItemSolver () override = default;

//    void
//    Compute (double time, bool view = true, bool force = false) override;


//    ItemSolver<t_type>&
//    SetEmbMesh(Mesh *mesh) override;

//    ItemSolver<t_type>&
//    SetTagToApply(PHYS tag) override;

//    ItemSolver<t_type>&
//    SetDerT(std::size_t dert) override;

//    ItemSolver<t_type>&
//    SetVariableOverTime(bool variableOverTime) override;

//    ItemSolver<t_type>&
//    SetFunction(std::function<double (Point, double)> fun) override;

//    template<ITEM_SOLVER_TYPE t_type_solver>
//    friend
//    std::ostream&
//    operator<< (std::ostream& out, const ItemSolver<t_type_solver>& item);
//};

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>::ItemSolver () :
//    ItemSolverBase ()
//{
//    m_type = t_type;
//}

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>::ItemSolver (const ItemSolver<t_type>& tocopy) :
//    ItemSolverBase (tocopy)
//{
//    m_type = t_type;
//}

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>&
//ItemSolver<t_type>::SetEmbMesh (Mesh* mesh)
//{
//    this->m_emb = mesh;
//    return *this;
//}

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>&
//ItemSolver<t_type>::SetTagToApply (PHYS tag)
//{
//    this->m_tagToApply = tag;
//    return *this;
//}

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>&
//ItemSolver<t_type>::SetFunction (std::function<double(Point, double)> fun)
//{
//    this->m_fun = fun;
//    return *this;
//}

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>&
//ItemSolver<t_type>::SetDerT (std::size_t dert)
//{
//    this->m_dert = dert;
//    return *this;
//}

//template <ITEM_SOLVER_TYPE t_type>
//ItemSolver<t_type>&
//ItemSolver<t_type>::SetVariableOverTime (bool variableOverTime)
//{
//    this->m_variableOverTime = variableOverTime;
//    return *this;
//}


//template <ITEM_SOLVER_TYPE t_type>
//std::ostream&
//operator<< (std::ostream& out, const ItemSolver<t_type>& item)
//{
//    out << "Item : "    << COLOR_GREEN  << std::setw (20) << item.m_type                                                   << FLUSHLINE;
//    out << " [der = "   << COLOR_BLUE   << item.m_dert                                                   << FLUSHLINE;
//    out << ", tag = "   << COLOR_BLUE   << to_string(item.m_tagToApply)                                             << FLUSHLINE;
//    out << ", pen = "   << COLOR_RED    << *item.m_coeff_pen << "/" << (*item.m_mesh)->GetPrescribedSize () << FLUSHLINE;
//    out << ", var = "   << COLOR_BLUE   << std::boolalpha << item.m_variableOverTime << std::noboolalpha << COLOR_DEFAULT << "]" << std::flush;

//    return out;
//}

//#define EXPLICIT_INST(X)                                                                            \
//    template<> void ItemSolver<ITEM_SOLVER_TYPE::X>::Compute (double time, bool view, bool force)

//EXPLICIT_INST (EMPTY);
//EXPLICIT_INST (PHI_PHI);
//EXPLICIT_INST (DAMPING_PHI_PHI);
//EXPLICIT_INST (GRAD_GRAD);
//EXPLICIT_INST (SECOND_MEMBER);
//EXPLICIT_INST (DIRICHLET_BOUNDARY);
//EXPLICIT_INST (NEUMANN_BOUNDARY);
//EXPLICIT_INST (DIRICHLET_BOUNDARY_SBM);
//EXPLICIT_INST (NEUMANN_BOUNDARY_SBM);

//#undef EXPLICIT_INST

//#endif // ITEMSOLVER_H
