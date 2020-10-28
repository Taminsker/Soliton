#include "superitem.h"

#include <Solver/SuperSolver/SuperSolver/supersolver.h>

double NULL_FUNC_NAME (Point p, double t)
{
    VOID_USE(p); VOID_USE(t);
    return 0.;
}

//SparseMatrix                            *m_mat;
//PlainVector                             *m_sec;
//std::vector<int>                        *m_appCells, *m_appEdges;



SuperItem::SuperItem (ITEM_T type) :
    m_type (type),
    m_mesh (nullptr),
    m_target (nullptr),
    m_solver (nullptr),
    m_mat(new SparseMatrix()),
    m_sec(new PlainVector ()),
    m_appCells (new std::vector<int>()),
    m_appEdges (new std::vector<int>()),
    m_appPoints (new std::vector<int>()),
    m_name (""),
    m_tag2Apply (PHYS::NONE),
    m_var0vrT (true),
    m_derT (0),
    m_fun (NULL_FUNC_NAME)
{}


SuperItem::SuperItem (const SuperItem& tocopy) :
    m_type (tocopy.m_type),
    m_mesh (nullptr),
    m_target (tocopy.m_target),
    m_solver (nullptr),
    m_mat (new SparseMatrix ()),
    m_sec (new PlainVector ()),
    m_appCells (new std::vector<int>()),
    m_appEdges (new std::vector<int>()),
    m_appPoints (new std::vector<int>()),
    m_name (""),
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
    delete m_appPoints;
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
    m_powpen                = &m_solver->m_pow_pen;
    m_time                  = &m_solver->m_time;

    std::size_t numCells    = static_cast<std::size_t>((*m_mesh)->GetNumberOfCells ());
    std::size_t numEdges    = static_cast<std::size_t>((*m_mesh)->GetNumberOfEdges ());
    m_appCells->resize (numCells, 0);
    m_appEdges->resize (numEdges, 0);

    return;
}

std::ostream& operator<< (std::ostream& out, const SuperItem& item)
{
    out << "Item : "    << COLOR_GREEN  << std::setw (22) << to_string(item.m_type) << FLUSHLINE;
    out << " [der = "   << COLOR_BLUE   << item.m_derT << FLUSHLINE;
    out << ", tags = "   << COLOR_BLUE   << std::setw (7) <<  to_string(item.m_tag2Apply) << COLOR_DEFAULT << "[" << COLOR_YELLOW << to_string(item.m_solver->GetInterTag ()) << COLOR_DEFAULT << "]" << FLUSHLINE;
    out << ", pen = "   << COLOR_RED    << *item.m_coeffpen << "/" << (*item.m_mesh)->GetPrescribedSize () << "^" << *item.m_powpen << FLUSHLINE;
    out << ", var = "   << COLOR_BLUE   << std::boolalpha << item.m_var0vrT << std::noboolalpha << COLOR_DEFAULT << "]" << std::flush;

    return out;
}
