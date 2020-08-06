#ifndef EDGE_H
#define EDGE_H

#include "../Cell/cell.h"

class Edge : public Cell
{
public:
    Edge ();
    ~Edge ();

    std::vector <Cell*>* GetCellsList ();
    int GetNumberOfCells () const;
    void AddCell (Cell* c);
    void RemoveCell (Cell* c);
    void RemoveCell (int idx);
    Cell* GetCell (int idx);
    void DetachFromAll ();

protected:
    std::vector <Cell*> m_celllist;

    std::vector <Edge *>* GetEdges () = delete;
    int GetNumberOfEdges () const = delete;
};

#endif // EDGE_H
