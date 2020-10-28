#include "mesh.h"
#include "../Cell/cell.h"
#include "../Edge/edge.h"
#include "../Point/point.h"
#include "../HetContainer/hetcontainer.h"

#include <ProgDef/proddef.h>
#include <Enums/enums.h>
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


void Mesh::Print () const
{
    BEGIN << "Mesh statistics." << ENDLINE;

    int count = 0;

    INFOS << COLOR_BLUE << "STATS : " << ENDLINE;
    INFOS << "Name              : \t\"" << COLOR_BLUE << m_name << COLOR_DEFAULT << "\"" << ENDLINE;
    INFOS << "H-Space           : \t" << m_h << ENDLINE;
    INFOS << "Number of Points  : \t" << m_points.size () << ENDLINE;
    INFOS << "Number of Cells   : \t" << m_cells.size () << ENDLINE;
    INFOS << "Number of Edges   : \t" << m_edges.size () << ENDLINE;

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "CELLS TYPE :" << ENDLINE;

    for (VTK_CELL_TYPE type : EnumClass<VTK_CELL_TYPE>())
    {
        count = CountCellType (type);
        if (count > 0)
            INFOS << "Number of " << to_string(type) << " : \t" << count << ENDLINE;
    }

    INFOS << SEPARATOR << ENDLINE;
    INFOS << COLOR_BLUE << "EDGES TYPE :" << ENDLINE;

    for (VTK_CELL_TYPE type : EnumClass<VTK_CELL_TYPE>())
    {
        count = CountEdgeType (type);
        if (count > 0)
            INFOS << "Number of " << to_string(type) << " : \t" << count << ENDLINE;
    }

    INFOS << SEPARATOR << ENDLINE;
    HetInt::Array* tagvec = GetPointsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec != nullptr)
    {
        INFOS << COLOR_BLUE << "POINTS TAG :" << ENDLINE;

        for (PHYS tag : EnumClass<PHYS>())
        {
            count = 0;
            for (auto tagid : tagvec->vec)
                if (static_cast<PHYS>(tagid) == tag)
                    count++;

            if (count > 0)
                INFOS << "Number of tag " << to_string(tag) << " : \t" << count << ENDLINE;
        }

        INFOS << SEPARATOR << ENDLINE;
    }


    tagvec = GetCellsData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec != nullptr)
    {
        INFOS << COLOR_BLUE << "CELLS TAG :" << ENDLINE;

        for (PHYS tag : EnumClass<PHYS>())
        {
            count = 0;
            for (auto tagid : tagvec->vec)
                if (static_cast<PHYS>(tagid) == tag)
                    count++;

            if (count > 0)
                INFOS << "Number of tag " << to_string(tag) << " : \t" << count << ENDLINE;
        }

        INFOS << SEPARATOR << ENDLINE;
    }

    tagvec = GetEdgesData ()->GetIntArrays ()->Get (NAME_TAG_PHYSICAL);

    if (tagvec != nullptr)
    {
        INFOS << COLOR_BLUE << "EDGES TAG :" << ENDLINE;

        for (PHYS tag : EnumClass<PHYS>())
        {
            count = 0;
            for (auto tagid : tagvec->vec)
                if (static_cast<PHYS>(tagid) == tag)
                    count++;

            if (count > 0)
                INFOS << "Number of tag " << to_string(tag) << " : \t" << count << ENDLINE;
        }

        INFOS << SEPARATOR << ENDLINE;
    }

    if (m_cellsdata->GetNumberOfArrays () != 0)
    {
        INFOS << COLOR_BLUE << "CELLS DATA :" << ENDLINE;
        m_cellsdata->Print ();
        INFOS << SEPARATOR << ENDLINE;
    }

    if (m_edgesdata->GetNumberOfArrays () != 0)
    {
        INFOS << COLOR_BLUE << "EDGES DATA :" << ENDLINE;
        m_edgesdata->Print ();
        INFOS << SEPARATOR << ENDLINE;
    }

    if (m_pointsdata->GetNumberOfArrays () != 0)
    {
        INFOS << COLOR_BLUE << "POINTS DATA :" << ENDLINE;
        m_pointsdata->Print ();
    }

    ENDFUN;
    return;
}
