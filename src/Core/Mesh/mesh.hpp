#ifndef SRC_CORE_MESH_MESH_HPP
#define SRC_CORE_MESH_MESH_HPP

#include <string>
#include <vector>

#include "../../Enums/enums.hpp"
#include "../../solitonheader.hpp"
#include "../Cell/cell.hpp"
#include "../Edge/edge.hpp"
#include "../HetContainer/hetcontainer.hpp"
#include "../Point/point.hpp"

struct InputDatStruct;

typedef DataContainer<Cell *>  CellsData;
typedef DataContainer<Edge *>  EdgesData;
typedef DataContainer<Point *> PointsData;

class Mesh
{
public:
    Mesh ();
    ~Mesh ();

    SOLITON_INLINE
    void
    AddPoint (Point * p)
    {
        m_points.push_back (p);
        return;
    }

    SOLITON_INLINE
    Point *
    GetPoint (ul_t i) const
    {
        return m_points [i];
    }

    SOLITON_INLINE
    void
    RemovePoint (ul_t i)
    {
        if (i < m_points.size ())
        {
            delete m_points [i];
            m_points.erase (m_points.begin () + i);
            m_pointsdata->Delete (i);
        }
        return;
    }

    SOLITON_INLINE
    void
    AddCell (Cell * c)
    {
        m_cells.push_back (c);
        return;
    }

    SOLITON_INLINE
    Cell *
    GetCell (ul_t i) const
    {
        return m_cells [i];
    }

    SOLITON_INLINE
    void
    RemoveCell (ul_t i)
    {
        if (i < m_cells.size ())
        {
            delete m_cells [i];
            m_cells.erase (m_cells.begin () + i);
            m_cellsdata->Delete (i);
        }

        return;
    }

    SOLITON_INLINE
    void
    AddEdge (Edge * c)
    {
        m_edges.push_back (c);
        return;
    }

    SOLITON_INLINE
    Edge *
    GetEdge (ul_t i) const
    {
        return m_edges [i];
    }

    SOLITON_INLINE
    void
    RemoveEdge (ul_t i)
    {
        if (i < m_edges.size ())
        {
            delete m_edges [i];
            m_edges.erase (m_edges.begin () + i);
            m_edgesdata->Delete (i);
        }

        return;
    }

    SOLITON_INLINE
    PointsData *
    GetPointsData () const
    {
        return m_pointsdata;
    }

    SOLITON_INLINE
    CellsData *
    GetCellsData () const
    {
        return m_cellsdata;
    }

    SOLITON_INLINE
    EdgesData *
    GetEdgesData () const
    {
        return m_edgesdata;
    }

    SOLITON_INLINE
    int
    GetNumberOfPoints () const
    {
        return int (m_points.size ());
    }

    SOLITON_INLINE
    int
    GetNumberOfCells () const
    {
        return int (m_cells.size ());
    }

    SOLITON_INLINE
    int
    GetNumberOfEdges () const
    {
        return int (m_edges.size ());
    }

    SOLITON_INLINE
    int
    GetNumberOfInfosCells () const
    {
        int sum = 0;
        for (ul_t i = 0; i < m_cells.size (); ++i)
            sum += m_cells [i]->GetNumberOfInfos ();
        return sum;
    }

    SOLITON_INLINE
    int
    GetNumberOfInfosEdges () const
    {
        int sum = 0;
        for (ul_t i = 0; i < m_edges.size (); ++i)
            sum += m_edges [i]->GetNumberOfInfos ();
        return sum;
    }

    void Print () const;

    SOLITON_INLINE
    void
    SetName (std::string name)
    {
        m_name = name;
        return;
    }

    SOLITON_INLINE
    std::string
    GetName () const
    {
        return m_name;
    }

    SOLITON_INLINE
    void
    SetPrescribedSize (real_t h)
    {
        m_h = h;
        return;
    }

    SOLITON_INLINE
    real_t
    GetPrescribedSize ()
    {
        return m_h;
    }

private:
    PointsData * m_pointsdata;
    CellsData *  m_cellsdata;
    EdgesData *  m_edgesdata;

    std::vector<Point *> m_points;
    std::vector<Cell *>  m_cells;
    std::vector<Edge *>  m_edges;

    std::string m_name;
    real_t      m_h;

    SOLITON_INLINE
    int
    CountCellType (VTK_CELL_TYPE type) const
    {
        int count = 0;
        for (auto c : m_cells)
            if (c->GetTypeVTK () == type)
                count++;
        return count;
    }

    SOLITON_INLINE
    int
    CountEdgeType (VTK_CELL_TYPE type) const
    {
        int count = 0;
        for (auto c : m_edges)
            if (c->GetTypeVTK () == type)
                count++;
        return count;
    }
};

#endif /* SRC_CORE_MESH_MESH_HPP */
