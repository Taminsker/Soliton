#include "cell.h"
#include "../Point/point.h"
#include "../Edge/edge.h"

#include <Enums/enums.h>

Cell::Cell() :
    m_globalIndex (NONE_ID_SELECTED), // Il n'a pas encore d'indice global
    m_centroid (new Point ()),
    m_points({}),
    m_edges ({})
{}

Cell::Cell (const Cell& c) :
    m_globalIndex (NONE_ID_SELECTED),
    m_cat_cell (c.m_cat_cell),
    m_celltag (c.m_celltag),
    m_celltype (c.m_celltype),
    m_centroid (new Point (*c.m_centroid)),
    m_points ({}),
    m_edges ({})
{
    for (Point* p : c.m_points)
        m_points.push_back (p);

    for (Edge* e : c.m_edges)
        m_edges.push_back (e);


    UpdateCentroid ();
}

Cell& Cell::operator= (std::initializer_list<Point *> ilist)
{
    m_points.clear ();

    for (auto iter = ilist.begin (); iter != ilist.end (); iter++)
        m_points.push_back (*iter);
    m_globalIndex = NONE_ID_SELECTED;

    UpdateCentroid ();

    return *this;
}

Cell::~Cell ()
{
    delete m_centroid;
    m_points.clear ();
    m_edges.clear ();
}


Cell* Cell::AddPoint (Point *p, bool l)
{
    for (ul_t i = 0; i < m_points.size (); ++i)
        if (m_points.at (i) == p)
            return this;

    m_points.push_back (p);
    if (l)
        p->LinkToCell (this);

    UpdateCentroid ();
    return this;
}

Cell* Cell::AddPoints (std::vector<Point*> plist)
{
    for (Point* p : plist)
    {
        for (ul_t i = 0; i < m_points.size (); ++i)
            if (m_points.at (i) == p)
                return this;

        m_points.push_back (p);
        p->LinkToCell (this);
    }

    UpdateCentroid ();
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

    UpdateCentroid ();
    return this;
}

Cell* Cell::AddEdge (Edge *e)
{
    for (ul_t i = 0; i < m_edges.size (); ++i)
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

void Cell::UpdateCentroid ()
{
    *m_centroid = {0, 0, 0};

    for (Point *p : m_points)
        *m_centroid += *p;

    *m_centroid = *m_centroid / static_cast<real_t> (m_points.size ());

    return;
}

std::ostream& operator<< (std::ostream &out, const Cell &c)
{
    out << c.m_points.size () << " ";
    for (ul_t i = 0; i < c.m_points.size (); ++i)
    {
        out << c.m_points.at (i)->GetGlobalIndex () << " ";
    }

    out << std::endl;

    return out;
}


