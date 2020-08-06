#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <vector>
#include <initializer_list>

#include "vtkCellType.h"

class Point;
class Edge;

typedef unsigned int CELLTAG;
typedef unsigned int CELLTYPE;
typedef enum {IMACELL, IMANEDGE} SPECIALCELL;

class Cell
{
public:
    Cell ();
    Cell (const Cell& c);
    ~Cell ();

    Cell& operator= (std::initializer_list<Point *> ilist);

    Cell* SetGlobalIndex (int index);
    int GetGlobalIndex () const;

    Cell* AddPoint (Point * p, bool link = true);
    Cell* RemovePoint (Point * p);

    Cell* AddEdge (Edge* e);
    Cell* RemoveEdge (Edge* e);

    void SetTag (CELLTAG tag);
    CELLTAG GetTag () const;

    void SetType (CELLTYPE type);

    CELLTYPE GetTypeGMSH () const;
    VTKCellType GetTypeVTK () const;

    int GetNumberOfInfos () const;

    std::vector <Point *>* GetPoints ();
    int GetNumberOfPoints () const;

    std::vector <Edge *>* GetEdges ();
    int GetNumberOfEdges () const;

    SPECIALCELL GetSpecial ();
    friend std::ostream & operator<< (std::ostream &out, const Cell &c);

protected:
    int m_globalIndex;
    SPECIALCELL m_specialcell = IMACELL;
    CELLTAG m_celltag;
    CELLTYPE m_celltype;
    std::vector <Point *> m_points;
    std::vector <Edge*> m_edges;
};

std::ostream & operator<< (std::ostream &out, const Cell &c);


#endif // CELL_H
