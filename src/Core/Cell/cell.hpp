#ifndef SRC_CORE_CELL_CELL_HPP
#define SRC_CORE_CELL_CELL_HPP

#include <fstream>
#include <initializer_list>
#include <vector>

#include "../../Enums/enums.hpp"
#include "../../solitonheader.hpp"

class Point;
class Edge;

typedef unsigned int CELLTAG;

class Cell
{
public:
    Cell ();
    Cell (const Cell & c);
    ~Cell ();

    Cell & operator= (std::initializer_list<Point *> ilist);

    Cell * AddPoint (Point * p, bool link = true);
    Cell * AddPoints (std::vector<Point *> plist);
    Cell * RemovePoint (Point * p);
    Cell * AddEdge (Edge * e);
    Cell * RemoveEdge (Edge * e);

    SOLITON_INLINE
    Cell *
    SetGlobalIndex (int index)
    {
        m_globalIndex = index;
        return this;
    }

    SOLITON_INLINE
    int
    GetGlobalIndex () const
    {
        return m_globalIndex;
    }

    SOLITON_INLINE
    void
    SetTag (CELLTAG tag)
    {
        m_celltag = tag;
        return;
    }

    SOLITON_INLINE
    CELLTAG
    GetTag () const
    {
        return m_celltag;
    }

    SOLITON_INLINE
    void
    SetType (GMSH_CELL_TYPE type)
    {
        m_celltype = type;
        return;
    }

    SOLITON_INLINE
    GMSH_CELL_TYPE
    GetTypeGMSH () const
    {
        return m_celltype;
    }

    SOLITON_INLINE
    VTK_CELL_TYPE
    GetTypeVTK () const
    {
        return Convert<GMSH_CELL_TYPE, VTK_CELL_TYPE> (m_celltype);
    }

    SOLITON_INLINE
    int
    GetNumberOfInfos () const
    {
        return 1 + int (m_points.size ());  // +1 pour le label de comptage de labels
    }

    SOLITON_INLINE
    std::vector<Point *> *
    GetPoints ()
    {
        return &m_points;
    }

    SOLITON_INLINE
    int
    GetNumberOfPoints () const
    {
        return static_cast<int> (m_points.size ());
    }

    SOLITON_INLINE
    std::vector<Edge *> *
    GetEdges ()
    {
        return &m_edges;
    }

    SOLITON_INLINE
    int
    GetNumberOfEdges () const
    {
        return static_cast<int> (m_edges.size ());
    }

    SOLITON_INLINE
    CAT_CELL_EDGE
    GetCat ()
    {
        return m_cat_cell;
    }

    SOLITON_INLINE
    Point *
    GetCentroid ()
    {
        return m_centroid;
    }

    SOLITON_INLINE
    void
    ForceUpdateCentroid ()
    {
        return UpdateCentroid ();
    }

    friend std::ostream & operator<< (std::ostream & out, const Cell & c);

protected:
    int                  m_globalIndex;
    CAT_CELL_EDGE        m_cat_cell = CAT_CELL_EDGE::CELL;
    CELLTAG              m_celltag;
    GMSH_CELL_TYPE       m_celltype;
    Point *              m_centroid;
    std::vector<Point *> m_points;
    std::vector<Edge *>  m_edges;

    void UpdateCentroid ();
};

std::ostream & operator<< (std::ostream & out, const Cell & c);

#endif /* SRC_CORE_CELL_CELL_HPP */
