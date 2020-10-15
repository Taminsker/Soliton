#include "edge.h"
#include "../Point/point.h"

Edge::Edge () :
    Cell(),
    m_celllist ({})
{
    this->m_cat_cell = CAT_CELL_EDGE::EDGE;
}

Edge::Edge (const Edge& tocopy) :
    Cell (tocopy),
    m_celllist ({})
{
    for (Cell* cell : tocopy.m_celllist)
        m_celllist.push_back (cell);

    this->m_cat_cell = CAT_CELL_EDGE::EDGE;
}

Edge::~Edge ()
{
    for (Cell* cell : m_celllist)
        cell->RemoveEdge (this);

    m_celllist.clear ();
}

void Edge::DetachFromAll ()
{
    for (Cell* c : m_celllist)
        c->RemoveEdge (this);

    for (Point* p : m_points)
        p->UnlinkToCell (this);

    return;
}
