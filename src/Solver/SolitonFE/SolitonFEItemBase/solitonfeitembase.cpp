#include "solitonfeitembase.hpp"

#include "../SolitonFEContainer/solitonfecontainer.hpp"

real_t NULL_FUNC_NAME (Point, real_t)
{
    return 0.;
}

SolitonFEItemBase::SolitonFEItemBase (SolitonFEContainer * parent, PHYS tag) : m_c (parent),
                                                                               m_id_target (0x0),
                                                                               m_mat (new SparseMatrix ()),
                                                                               m_sec (new DenseVector ()),
                                                                               m_appCells (new std::vector<int> ()),
                                                                               m_appEdges (new std::vector<int> ()),
                                                                               m_appPoints (new std::vector<int> ()),
                                                                               m_name (""),
                                                                               m_tag2Apply (tag),
                                                                               m_var0vrT (false),
                                                                               m_derT (0),
                                                                               m_fun (&NULL_SOLITON_FUNCTOR)
{
    ul_t numCells = static_cast<ul_t> (m_c->m_store->GetMainMesh ()->GetNumberOfCells ());
    ul_t numEdges = static_cast<ul_t> (m_c->m_store->GetMainMesh ()->GetNumberOfEdges ());
    m_appCells->resize (numCells, 0);
    m_appEdges->resize (numEdges, 0);
}

SolitonFEItemBase::SolitonFEItemBase (const SolitonFEItemBase & tocopy) : m_c (tocopy.m_c),
                                                                          m_id_target (tocopy.m_id_target),
                                                                          m_mat (new SparseMatrix ()),
                                                                          m_sec (new DenseVector ()),
                                                                          m_appCells (new std::vector<int> ()),
                                                                          m_appEdges (new std::vector<int> ()),
                                                                          m_appPoints (new std::vector<int> ()),
                                                                          m_name (tocopy.m_name),
                                                                          m_tag2Apply (tocopy.m_tag2Apply),
                                                                          m_var0vrT (tocopy.m_var0vrT),
                                                                          m_derT (tocopy.m_derT),
                                                                          m_fun (tocopy.m_fun)
{
    ul_t numCells = static_cast<ul_t> (m_c->m_store->GetMainMesh ()->GetNumberOfCells ());
    ul_t numEdges = static_cast<ul_t> (m_c->m_store->GetMainMesh ()->GetNumberOfEdges ());
    m_appCells->resize (numCells, 0);
    m_appEdges->resize (numEdges, 0);
}

SolitonFEItemBase::~SolitonFEItemBase ()
{
    delete m_mat;
    delete m_sec;
    delete m_appCells;
    delete m_appEdges;
    delete m_appPoints;
}

void
SolitonFEItemBase::Print (std::ostream & out)
{
    out << "Item : " << COLOR_GREEN << std::setw (22) << to_string (GetType ()) << FLUSHLINE;
    out << " [der = " << COLOR_BLUE << m_derT << FLUSHLINE;
    out << ", tags = " << COLOR_BLUE << std::setw (7) << to_string (m_tag2Apply) << COLOR_DEFAULT << "[" << COLOR_YELLOW << to_string (m_c->GetInterTag ()) << COLOR_DEFAULT << "]" << FLUSHLINE;
    out << ", pen = " << COLOR_RED << m_c->GetCoeffPen () << "/" << m_c->m_store->GetMainMesh ()->GetPrescribedSize () << "^" << m_c->GetPowPen () << FLUSHLINE;
    out << ", var = " << COLOR_BLUE << std::boolalpha << m_var0vrT << std::noboolalpha << COLOR_DEFAULT << "]" << std::flush;

    return;
}

void
SolitonFEItemBase::SetSparseMatrix (std::vector<Triplet> * listTriplets)
{
    int numPoints = m_c->m_store->GetMainMesh ()->GetNumberOfPoints ();
    m_mat->resize (numPoints, numPoints);
    m_mat->setFromTriplets (listTriplets->begin (), listTriplets->end ());
    return;
}

void
SolitonFEItemBase::SetSecMember (std::vector<Duo> * listDuos)
{
    int numPoints = m_c->m_store->GetMainMesh ()->GetNumberOfPoints ();

    m_sec->setZero (numPoints);

    for (Duo duo : *listDuos)
        m_sec->coeffRef (duo.idx ()) += duo.value ();

    return;
}
