#include "febase.h"

FEBase::FEBase () :
    m_type_vtk (VTK_CELL_TYPE::VTK_EMPTY_CELL),
    m_volref (0.0),
    m_order (0),
    m_dim(0),
    m_type_physical(PHYSICAL_CELL_TYPE::EMPTY),
    m_coor ({}),
    m_edges ({}),
    m_edges_normals ({}),
    m_edges_tangents ({}),
    m_phi ({}),
    m_grad_phi ({}),
    m_target (nullptr),
    m_cat (FE_CLASS_TYPE::EMPTY)

{}

FEBase::~FEBase ()
{
    m_phi.clear ();
    m_grad_phi.clear ();

    for (Point* p : m_coor)
        delete p;
    m_coor.clear ();

    for (Edge* e : m_edges)
        delete e;
    m_edges.clear ();

    for (Point* p : m_edges_normals)
        delete p;
    m_edges_normals.clear ();

    for (Point* p : m_edges_tangents)
        delete  p;
    m_edges_tangents.clear ();
}

void FEBase::CastForCell (Cell* target)
{
    this->m_target = target;

    return;
}

Cell* FEBase::GetTarget ()
{
    return m_target;
}


VTK_CELL_TYPE FEBase::GetCellType ()
{
    return m_type_vtk;
}

std::size_t FEBase::GetNumberOfPoints ()
{
    return m_coor.size ();
}

std::size_t FEBase::GetNumberOfEdges ()
{
    return m_edges.size ();
}

double FEBase::GetMesureOfRef ()
{
    return m_volref;
}

std::size_t FEBase::GetOrderOfElement ()
{
    return m_order;
}

std::size_t FEBase::GetDimension ()
{
    return m_dim;
}

PHYSICAL_CELL_TYPE FEBase::GetTypePhysicalBase ()
{
    return m_type_physical;
}

Point* FEBase::GetPoint (std::size_t id)
{
    return m_coor.at (id);
}

Edge* FEBase::GetEdge (std::size_t id)
{
    return m_edges.at (id);
}

Point* FEBase::GetNormalToEdge (std::size_t id)
{
    return m_edges_normals.at (id);
}

Point* FEBase::GetTangentToEdge (std::size_t id)
{
    return m_edges_tangents.at (id);
}

LambdaOnPoint<double>& FEBase::GetPhi (std::size_t id)
{
    return m_phi.at (id);
}

LambdaOnPoint<Point>& FEBase::GetGradPhi (std::size_t id)
{
    return m_grad_phi.at (id);
}
