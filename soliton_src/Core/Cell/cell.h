#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <vector>
#include <initializer_list>

#include "celltype.h"

class Point;
class Edge;

typedef unsigned int CELLTAG;
enum class CAT_CELL_EDGE
{
    CELL    = 0,
    EDGE    = 1,
    FIRST = CELL,
    LAST = EDGE
};

class Cell
{
public:
    Cell ();
    Cell (const Cell& c);
    ~Cell ();

    Cell&                   operator= (std::initializer_list<Point *> ilist);

    Cell*                   SetGlobalIndex (int index);
    int                     GetGlobalIndex () const;
    Cell*                   AddPoint (Point * p, bool link = true);
    Cell*                   AddPoints (std::vector<Point*> plist);
    Cell*                   RemovePoint (Point * p);
    Cell*                   AddEdge (Edge* e);
    Cell*                   RemoveEdge (Edge* e);
    void                    SetTag (CELLTAG tag);
    CELLTAG                 GetTag () const;
    void                    SetType (GMSH_CELL_TYPE type);
    GMSH_CELL_TYPE          GetTypeGMSH () const;
    VTK_CELL_TYPE           GetTypeVTK () const;
    int                     GetNumberOfInfos () const;
    std::vector <Point *>*  GetPoints ();
    int                     GetNumberOfPoints () const;
    std::vector <Edge *>*   GetEdges ();
    int                     GetNumberOfEdges () const;
    CAT_CELL_EDGE           GetCat ();
    friend std::ostream&    operator<< (std::ostream &out, const Cell &c);

    Point* GetCentroid ();

protected:
    int                     m_globalIndex;
    CAT_CELL_EDGE           m_cat_cell = CAT_CELL_EDGE::CELL;
    CELLTAG                 m_celltag;
    GMSH_CELL_TYPE          m_celltype;
    Point*                  m_centroid;
    std::vector <Point*>    m_points;
    std::vector <Edge*>     m_edges;

    void                    UpdateCentroid ();
};

std::ostream& operator<< (std::ostream &out, const Cell &c);


#endif // CELL_H
