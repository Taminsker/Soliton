#include "edge.h"
#include "../Point/point.h"

Edge::Edge () :
    Cell(),
    m_celllist ({})
{
    this->m_cat_cell = CAT_CELL_EDGE::EDGE;
}

Edge::~Edge ()
{
    for (Cell* cell : m_celllist)
        cell->RemoveEdge (this);

    m_celllist.clear ();
}

std::vector<Cell*>* Edge::GetCellsList ()
{
    return &m_celllist;
}

int Edge::GetNumberOfCells () const
{
    return static_cast<int>(m_celllist.size ());
}

Cell* Edge::GetCell (int idx)
{
    if (idx >= 0 && idx < static_cast<int>(m_celllist.size ()))
        return m_celllist.at (std::size_t(idx));
    return nullptr;
}

void Edge::DetachFromAll ()
{
    for (Cell* c : m_celllist)
        c->RemoveEdge (this);

    for (Point* p : m_points)
        p->UnlinkToCell (this);

    return;
}

void Edge::AddCell (Cell *c)
{
    for (Cell* cell : m_celllist)
        if (c == cell)
            return;
    m_celllist.push_back (c);
    c->AddEdge (this);

    return;
}

void Edge::RemoveCell (Cell *c)
{
    auto it = m_celllist.begin ();
    while (it != m_celllist.end ())
    {
        if (*it == c)
        {
            (*it)->RemoveEdge (this);
            it = m_celllist.erase (it);
            return;
        }
        else
            it++;
    }
    return;
}


void Edge::RemoveCell (int idx)
{
    return RemoveCell (GetCell (idx));
}
