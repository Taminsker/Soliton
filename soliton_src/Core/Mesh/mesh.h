#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>

#include "../Cell/celltype.h"
#include "tagphysical.h"
#include <ProgDef/proddef.h>

struct InputDatStruct;
class Cell;
class Point;
class Edge;

template <typename T>
class DataContainer;

typedef DataContainer <Cell*>  CellsData;
typedef DataContainer <Edge*>  EdgesData;
typedef DataContainer <Point*> PointsData;

class Mesh
{
public:
    Mesh ();
    ~Mesh ();

    void            AddPoint (Point * p);
    Point*          GetPoint (int i) const;
    void            RemovePoint (int i);

    void            AddCell (Cell* c);
    Cell*           GetCell (int i) const;
    void            RemoveCell (int i);

    void            AddEdge (Edge* c);
    Edge*           GetEdge (int i) const;
    void            RemoveEdge (int i);

    PointsData*     GetPointsData () const;
    CellsData*      GetCellsData () const;
    EdgesData*      GetEdgesData () const;

    int             GetNumberOfPoints () const;
    int             GetNumberOfCells () const;
    int             GetNumberOfEdges () const;

    int             GetNumberOfInfosCells () const;
    int             GetNumberOfInfosEdges () const;

    void            Print () const;

    void            SetName (std::string name);
    std::string     GetName () const;

    void            SetPrescribedSize(double h);
    double          GetPrescribedSize ();

private:
    PointsData*                     m_pointsdata;
    CellsData*                      m_cellsdata;
    EdgesData*                      m_edgesdata;

    std::vector <Point*>            m_points;
    std::vector <Cell*>             m_cells;
    std::vector <Edge*>             m_edges;

    std::string                     m_name;
    double                          m_h;


    int             CountCellType (VTK_CELL_TYPE type) const;
    int             CountEdgeType (VTK_CELL_TYPE type) const;

};

#endif // MESH_H
