#include "superitem.h"

#include <Solver/SuperSolver/SuperSolver/supersolver.h>

SuperItem::SuperItem (ITEM_T type) :
    m_type (type),
    m_mat(new SparseMatrix()),
    m_sec(new PlainVector ()),
    m_tag2Apply (PHYS::NONE),
    m_var0vrT (true),
    m_derT (0),
    m_fun (NULL_FUNC_NAME)
{}


SuperItem::SuperItem (const SuperItem& tocopy) :
    m_type (tocopy.m_type),
    m_mat (new SparseMatrix ()),
    m_sec (new PlainVector ()),
    m_tag2Apply (tocopy.m_tag2Apply),
    m_var0vrT (tocopy.m_var0vrT),
    m_derT (tocopy.m_derT),
    m_fun (tocopy.m_fun)
{
    this->SetSuperSolver (tocopy.m_solver);
}


SuperItem::~SuperItem ()
{
    delete m_mat;
    delete m_sec;
    delete m_appCells;
    delete m_appEdges;
}

void SuperItem::SetSuperSolver (SuperSolver *solver)
{
    if (solver == nullptr)
        return;

    m_mesh                  = &solver->m_store->mesh;
    m_solver                = solver;
    m_collect               = &m_solver->m_collect_contributions;
    m_view                  = &m_solver->m_view;
    m_force_compute         = &m_solver->m_force_compute;
    m_coeffpen              = &m_solver->m_coeff_pen;
    m_time                  = &m_solver->m_time;

    std::size_t numCells    = static_cast<std::size_t>(m_solver->m_store->mesh->GetNumberOfCells ());
    std::size_t numEdges    = static_cast<std::size_t>(m_solver->m_store->mesh->GetNumberOfEdges ());
    m_appCells              = new std::vector<int>(numCells, 0);
    m_appEdges              = new std::vector<int>(numEdges, 0);

    return;
}

std::ostream& operator<< (std::ostream& out, const SuperItem& item)
{
    out << "Item : "    << COLOR_GREEN  << to_string(item.m_type) << FLUSHLINE;
    out << " [der = "   << COLOR_BLUE   << item.m_derT << FLUSHLINE;
    out << ", tag = "   << COLOR_BLUE   << to_string(item.m_tag2Apply) << FLUSHLINE;
    out << ", pen = "   << COLOR_RED    << *item.m_coeffpen << "/" << (*item.m_mesh)->GetPrescribedSize () << FLUSHLINE;
    out << ", var = "   << COLOR_BLUE   << std::boolalpha << item.m_var0vrT << std::noboolalpha << COLOR_DEFAULT << "]" << std::flush;

    return out;
}
