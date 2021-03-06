#ifndef SRC_CORE_EDGE_EDGE_HPP
#define SRC_CORE_EDGE_EDGE_HPP

#include "../../solitonheader.hpp"
#include "../Cell/cell.hpp"

class Edge : public Cell
{
public:
    Edge ();
    Edge (const Edge & tocopy);
    ~Edge ();

    SOLITON_INLINE
    std::vector<Cell *> *
    GetCellsList ()
    {
        return &m_celllist;
    }

    SOLITON_INLINE
    int
    GetNumberOfCells () const
    {
        return static_cast<int> (m_celllist.size ());
    }

    SOLITON_INLINE
    void
    AddCell (Cell * c)
    {
        for (Cell * cell : m_celllist)
            if (c == cell)
                return;

        m_celllist.push_back (c);

        c->AddEdge (this);

        return;
    }

    SOLITON_INLINE
    void
    RemoveCell (Cell * c)
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

    SOLITON_INLINE
    void
    RemoveCell (int idx)
    {
        return RemoveCell (GetCell (idx));
    }

    SOLITON_INLINE
    Cell *
    GetCell (int idx)
    {
        if (idx >= 0 && idx < static_cast<int> (m_celllist.size ()))
            return m_celllist.at (ul_t (idx));
        return nullptr;
    }

    void DetachFromAll ();

protected:
    std::vector<Cell *> m_celllist;

    std::vector<Edge *> * GetEdges ()               = delete;
    int                   GetNumberOfEdges () const = delete;
};

#endif /* SRC_CORE_EDGE_EDGE_HPP */
