#include "cell.h"
#include "../Point/point.h"
#include "../Edge/edge.h"

#include "vtkCellType.h"

Cell::Cell() :
    m_globalIndex (NONEID), // Il n'a pas encore d'indice global
    m_points({}),
    m_edges ({})
{}

Cell::Cell (const Cell& c) :
    m_globalIndex (NONEID),
    m_points (c.m_points),
    m_edges (c.m_edges)
{}

Cell& Cell::operator= (std::initializer_list<Point *> ilist)
{
    m_points.clear ();

    for (auto iter = ilist.begin (); iter != ilist.end (); iter++)
        m_points.push_back (*iter);
    m_globalIndex = NONEID;

    return *this;
}

Cell::~Cell ()
{
    m_points.clear ();
    m_edges.clear ();
}


Cell* Cell::SetGlobalIndex (int index)
{
    m_globalIndex = index;
    return this;
}

int Cell::GetGlobalIndex () const
{
    return m_globalIndex;
}

Cell* Cell::AddPoint (Point *p, bool l)
{
    for (std::size_t i = 0; i < m_points.size (); ++i)
        if (m_points.at (i) == p)
            return this;

    m_points.push_back (p);
    if (l)
        p->LinkToCell (this);
    return this;
}

Cell* Cell::RemovePoint (Point *p)
{
    auto it = m_points.begin ();
    while (it != m_points.end ())
    {
        if (*it == p)
        {
            (*it)->UnlinkToCell (this);
            it = m_points.erase (it);
        }
        else
            it++;
    }

    return this;
}

Cell* Cell::AddEdge (Edge *e)
{
    for (std::size_t i = 0; i < m_edges.size (); ++i)
        if (m_edges.at (i) == e)
            return this;

    m_edges.push_back (e);
    return this;
}

Cell* Cell::RemoveEdge (Edge *e)
{
    auto it = m_edges.begin ();
    while (it != m_edges.end ())
    {
        if (*it == e)
        {
            it = m_edges.erase (it);
        }
        else
            it++;
    }

    return this;
}

void Cell::SetTag (CELLTAG tag)
{
    m_celltag = tag;
    return;
}

CELLTAG Cell::GetTag () const
{
    return m_celltag;
}

void Cell::SetType (CELLTYPE type)
{
    m_celltype = type;
    return;
}

CELLTYPE Cell::GetTypeGMSH () const
{
    return m_celltype;
}

VTKCellType Cell::GetTypeVTK () const
{
    switch (m_celltype) {
    default:
    case 0:
        return VTK_EMPTY_CELL;
    case 1: // 2-node line.
        return VTK_LINE;
    case 2: // 3-node triangle.
        return VTK_TRIANGLE;
    case 3: // 4-node quadrangle.
        return VTK_QUAD;
    case 4: // 4-node tetrahedron.
        return VTK_TETRA;
    case 5: // 8-node hexahedron.
        return VTK_HEXAHEDRON;
        //    case 6: // 6-node prism.
        //        return VTK_PYRAMID;
    case 7: // 5-node pyramid.
        return VTK_QUADRATIC_PYRAMID;
    case 8: // 3-node second order line
        return VTK_QUADRATIC_EDGE;
    case 9: // 6-node second order triangle
        return VTK_QUADRATIC_TRIANGLE;
    case 10: // 9-node second order quadrangle
        return VTK_QUADRATIC_QUAD;
    case 11: // 10-node second order tetrahedron
        return VTK_QUADRATIC_TETRA;
    case 12: // 27-node second order hexahedron
        return VTK_QUADRATIC_HEXAHEDRON;
        //    case 13: // 18-node second order prism
        //        return ;       //
    case 14: // 14-node second order pyramid
        return VTK_QUADRATIC_PYRAMID;
    case 15: // 1-node point.
        return VTK_VERTEX;
    case 16: // 8-node second order quadrangle
        return VTK_BIQUADRATIC_QUAD;
    case 17: // 20-node second order hexahedron
        return VTK_TRIQUADRATIC_HEXAHEDRON;
        //    case 18:       // 15-node second order prism
        //        return ;       //
        //    case 19:       // 13-node second order pyramid
        //        return ;       //
    }
}

int Cell::GetNumberOfInfos () const
{
    return 1 + int (m_points.size ()); // +1 pour le label de comptage de labels
}

std::vector <Point *>* Cell::GetPoints ()
{
    return &m_points;
}
std::vector <Edge *>* Cell::GetEdges ()
{
    return &m_edges;
}

int Cell::GetNumberOfPoints () const
{
    return int(m_points.size ());
}

int Cell::GetNumberOfEdges () const
{
    return int(m_edges.size ());
}

SPECIALCELL Cell::GetSpecial ()
{
    return m_specialcell;
}

std::ostream& operator<< (std::ostream &out, const Cell &c)
{
    out << c.m_points.size () << " ";
    for (std::size_t i = 0; i < c.m_points.size (); ++i)
    {
        out << c.m_points.at (i)->GetGlobalIndex () << " ";
    }

    out << std::endl;

    return out;
}
