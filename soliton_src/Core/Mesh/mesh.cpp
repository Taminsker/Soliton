#include "mesh.h"
#include "../Cell/cell.h"
#include "../Edge/edge.h"
#include "../Point/point.h"
#include "../HetContainer/hetcontainer.h"

#include <ProgDef/proddef.h>
#include <Core/EnumClass/enumclass.h>
#include <IO/parsers.h>

Mesh::Mesh () :
    m_pointsdata (new PointsData (&m_points)),
    m_cellsdata (new CellsData (&m_cells)),
    m_edgesdata (new EdgesData (&m_edges)),
    m_name ("no-name-selected"),
    m_h (1.)
{}

Mesh::~Mesh ()
{
    for (auto e : m_edges)
        delete e;
    for (auto p : m_points)
        delete p;
    for (auto c : m_cells)
        delete c;

    m_cells.clear();
    m_points.clear ();
    m_edges.clear ();

    delete m_cellsdata;
    delete m_edgesdata;
    delete m_pointsdata;
}

Point* Mesh::GetPoint (int index) const
{
    return m_points.at (std::size_t (index));
}

Cell* Mesh::GetCell (int index) const
{
    return m_cells.at (std::size_t (index));
}

Edge* Mesh::GetEdge (int index) const
{
    return m_edges.at (std::size_t (index));
}

PointsData* Mesh::GetPointsData () const
{
    return m_pointsdata;
}

CellsData* Mesh::GetCellsData () const
{
    return m_cellsdata;
}

EdgesData* Mesh::GetEdgesData () const
{
    return m_edgesdata;
}

void Mesh::AddPoint (Point * p)
{
    m_points.push_back (p);
    return;
}

void Mesh::AddCell (Cell* c)
{
    m_cells.push_back (c);
    return;
}

void Mesh::AddEdge (Edge* c)
{
    m_edges.push_back (c);
    return;
}

void Mesh::RemovePoint (int i)
{
    if (i >= 0 && i < int (m_points.size ()))
    {
        delete m_points.at (std::size_t (i));
        m_points.erase (m_points.begin () + i);
        m_pointsdata->Delete (i);
    }

    return;
}

void Mesh::RemoveCell (int i)
{
    if (i >= 0 && i < int (m_cells.size ()))
    {
        delete m_cells.at (std::size_t (i));
        m_cells.erase (m_cells.begin () + i);
        m_cellsdata->Delete (i);
    }

    return;
}

void Mesh::RemoveEdge (int i)
{
    if (i >= 0 && i < int (m_edges.size ()))
    {
        delete m_edges.at (std::size_t (i));
        m_edges.erase (m_edges.begin () + i);
        m_edgesdata->Delete (i);
    }

    return;
}

int Mesh::GetNumberOfPoints () const
{
    return int (m_points.size ());
}

int Mesh::GetNumberOfCells () const
{
    return int (m_cells.size ());
}

int Mesh::GetNumberOfEdges () const
{
    return int (m_edges.size ());
}

int Mesh::GetNumberOfInfosCells () const
{
    int sum = 0;
    for (std::size_t i = 0; i < m_cells.size (); ++i)
        sum += m_cells.at (i)->GetNumberOfInfos ();
    return sum;
}

int Mesh::GetNumberOfInfosEdges () const
{
    int sum = 0;
    for (std::size_t i = 0; i < m_edges.size (); ++i)
        sum += m_edges.at (i)->GetNumberOfInfos ();
    return sum;
}

void Mesh::SetName (std::string name)
{
    m_name = name;
    return;
}

std::string Mesh::GetName () const
{
    return m_name;
}

void Mesh::SetPrescribedSize(double h)
{
    m_h = h;
    return;
}

double Mesh::GetPrescribedSize ()
{
    return m_h;
}

void Mesh::Print () const
{
    HEADERFUN("Mesh::Print");
    BEGIN << "Mesh statistics." << ENDLINE;

    int count = 0;

    INFOS << COLOR_BLUE << "STATS : " << ENDLINE;
    INFOS << "Name :\t \"" << m_name << "\"" << ENDLINE;
    INFOS << "Number of Points : \t" << m_points.size () << ENDLINE;
    INFOS << "Number of Cells  : \t" << m_cells.size () << ENDLINE;
    INFOS << "Number of Edges  : \t" << m_edges.size () << ENDLINE;

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "CELLS TYPE :" << ENDLINE;

    for (VTK_CELL_TYPE type : EnumClass<VTK_CELL_TYPE>())
    {
        count = CountCellType (type);
        if (count > 0)
            INFOS << "Number of " << type << " : \t" << count << ENDLINE;
    }

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "EDGES TYPE :" << ENDLINE;

    for (VTK_CELL_TYPE type : EnumClass<VTK_CELL_TYPE>())
    {
        count = CountEdgeType (type);
        if (count > 0)
            INFOS << "Number of " << type << " : \t" << count << ENDLINE;
    }

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "POINTS TAG :" << ENDLINE;
    HetInt::Array* tagvec = GetPointsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec != nullptr)
    {
        for (TAG_PHYSICAL tag : EnumClass<TAG_PHYSICAL>())
        {
            count = 0;
            for (auto tagid : tagvec->vec)
                if (static_cast<TAG_PHYSICAL>(tagid) == tag)
                    count++;

            if (count > 0)
                INFOS << "Number of tag " << tag << " : \t" << count << ENDLINE;
        }
    }

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "CELLS TAG :" << ENDLINE;
    tagvec = GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec != nullptr)
    {
        for (TAG_PHYSICAL tag : EnumClass<TAG_PHYSICAL>())
        {
            count = 0;
            for (auto tagid : tagvec->vec)
                if (static_cast<TAG_PHYSICAL>(tagid) == tag)
                    count++;

            if (count > 0)
                INFOS << "Number of tag " << tag << " : \t" << count << ENDLINE;
        }
    }

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "EDGES TAG :" << ENDLINE;
    tagvec = GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec != nullptr)
    {
        for (TAG_PHYSICAL tag : EnumClass<TAG_PHYSICAL>())
        {
            count = 0;
            for (auto tagid : tagvec->vec)
                if (static_cast<TAG_PHYSICAL>(tagid) == tag)
                    count++;

            if (count > 0)
                INFOS << "Number of tag " << tag << " : \t" << count << ENDLINE;
        }
    }

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "CELLS DATA :" << ENDLINE;
    m_cellsdata->Print ();

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "EDGES DATA :" << ENDLINE;
    m_edgesdata->Print ();

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "POINTS DATA :" << ENDLINE;
    m_pointsdata->Print ();

    ENDFUN;
    return;
}

int Mesh::CountCellType (VTK_CELL_TYPE type) const
{
    int count = 0;
    for (auto c : m_cells)
        if (c->GetTypeVTK () == type)
            count++;
    return count;
}

int Mesh::CountEdgeType (VTK_CELL_TYPE type) const
{
    int count = 0;
    for (auto c : m_edges)
        if (c->GetTypeVTK () == type)
            count++;
    return count;
}

