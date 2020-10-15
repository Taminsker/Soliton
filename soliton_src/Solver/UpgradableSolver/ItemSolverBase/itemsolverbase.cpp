//#include "itemsolverbase.h"

//double NULL_FUNC_NAME (Point p, double t)
//{
//    VOID_USE(p);
//    VOID_USE(t);

//    return 0.;
//}

//ItemSolverBase::ItemSolverBase () :
//    m_mesh (nullptr),
//    m_coeff_pen (nullptr),
//    m_dfstore (new DFStore ()),
//    m_festore (new FEStore ()),
//    m_quadstore (new QuadStore ()),
//    m_triplets ({}),
//    m_duos ({}),
//    m_applicableCells ({}),
//    m_applicableEdges ({}),
//    m_emb (nullptr),
//    m_fun (nullptr),
//    m_tagToApply (PHYS::NONE),
//    m_dert (0),
//    m_type (ITEM_SOLVER_TYPE::EMPTY),
//    m_variableOverTime (false),
//    m_mat (),
//    m_second_member ()
//{}

//ItemSolverBase::ItemSolverBase (const ItemSolverBase& tocopy) :
//    m_mesh (tocopy.m_mesh),
//    m_coeff_pen (tocopy.m_coeff_pen),
//    m_dfstore (new DFStore ()),
//    m_festore (new FEStore ()),
//    m_quadstore (new QuadStore ()),
//    m_triplets ({}),
//    m_duos ({}),
//    m_applicableCells ({}),
//    m_applicableEdges ({}),
//    m_emb (tocopy.m_emb),
//    m_fun (tocopy.m_fun),
//    m_tagToApply (tocopy.m_tagToApply),
//    m_dert (tocopy.m_dert),
//    m_type (tocopy.m_type),
//    m_variableOverTime (tocopy.m_variableOverTime),
//    m_mat (),
//    m_second_member ()
//{}

//ItemSolverBase::~ItemSolverBase ()
//{
//    delete m_dfstore;
//    delete m_festore;
//    delete m_quadstore;
//}

//ItemSolverBase::operator SparseMatrix* ()
//{
//    return this->CastSparseMatrix ();
//}

//ItemSolverBase::operator PlainVector* ()
//{
//    return this->CastPlainVector ();
//}

//SparseMatrix*
//ItemSolverBase::CastSparseMatrix ()
//{
//    int numPoints = (*m_mesh)->GetNumberOfPoints ();
//    m_mat.resize (numPoints, numPoints);
//    m_mat.setFromTriplets (m_triplets.begin (), m_triplets.end ());

//    return &m_mat;
//}

//PlainVector*
//ItemSolverBase::CastPlainVector ()
//{
//    int numPoints = (*m_mesh)->GetNumberOfPoints ();
//    m_second_member.resize (numPoints);
//    m_second_member.setZero ();

//    for (Duo duo : m_duos)
//        m_second_member.coeffRef (long(duo.row ())) += duo.value ();

//    return &m_second_member;
//}

//std::vector<int>*
//ItemSolverBase::GetApplicableCellsVector ()
//{
//    return &m_applicableCells;
//}

//std::vector<int>*
//ItemSolverBase::GetApplicableEdgesVector ()
//{
//    return &m_applicableEdges;
//}

//void
//ItemSolverBase::ConfigureFromSolver (Mesh **mesh, double *coeffpen, bool *collect)
//{
//    m_mesh = mesh;
//    m_coeff_pen = coeffpen;
//    m_collect = collect;
//    return;
//}

//std::size_t
//ItemSolverBase::GetTemporalDerivationOrder ()
//{
//    return this->m_dert;
//}

//ITEM_SOLVER_TYPE
//ItemSolverBase::GetTypeOfItem ()
//{
//    return m_type;
//}

//bool
//ItemSolverBase::ItsVariableOverTime ()
//{
//    return this->m_variableOverTime;
//}

//PHYS
//ItemSolverBase::GetTagToApply ()
//{
//    return m_tagToApply;
//}


//std::ostream&
//operator<< (std::ostream& out, const ItemSolverBase& item)
//{
//    out << "Item : "    << COLOR_GREEN  << item.m_type                                                   << FLUSHLINE;
//    out << " [der = "   << COLOR_BLUE   << item.m_dert                                                   << FLUSHLINE;
//    out << ", tag = "   << COLOR_BLUE   << to_string(item.m_tagToApply)                                             << FLUSHLINE;
//    out << ", pen = "   << COLOR_RED    << *item.m_coeff_pen << "/" << (*item.m_mesh)->GetPrescribedSize () << FLUSHLINE;
//    out << ", var = "   << COLOR_BLUE   << std::boolalpha << item.m_variableOverTime << std::noboolalpha << COLOR_DEFAULT << "]" << std::flush;

//    return out;
//}

