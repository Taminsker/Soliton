#include "festore.h"

FEStore::FEStore() :
    m_emp0n0ddl (new Emp0N0DDL ()),
    m_ver1n1ddl (new Ver1N1DDL ()),
    m_lin2n2ddl (new Lin2N2DDL ()),
    m_lin3n3ddl (new Lin3N3DDL ()),
    m_tri3n3ddl (new Tri3N3DDL ()),
    m_tri6n6ddl (new Tri6N6DDL ()),
    m_quad4n4ddl (new Quad4N4DDL ()),
    m_quad8n8ddl (new Quad8N8DDL ()),
    m_quad9n9ddl (new Quad9N9DDL ()),
    m_tet4n4ddl (new Tet4N4DDL ())
{}

FEStore::~FEStore ()
{
    delete m_emp0n0ddl;
    delete m_ver1n1ddl;
    delete m_lin2n2ddl;
    delete m_lin3n3ddl;
    delete m_tri3n3ddl;
    delete m_tri6n6ddl;
    delete m_quad4n4ddl;
    delete m_quad8n8ddl;
    delete m_quad9n9ddl;
    delete m_tet4n4ddl;
}

FEBase* FEStore::GetElementFor (Cell* cell, FE_CLASS_TYPE typeOfEle)
{
    if (typeOfEle == FE_CLASS_TYPE::LAGRANGE)
    {
        switch (cell->GetTypeVTK ())
        {
        case VTK_CELL_TYPE::VTK_VERTEX:
            m_ver1n1ddl->CastForCell (cell);
            return m_ver1n1ddl;
        case VTK_CELL_TYPE::VTK_LINE:
            m_lin2n2ddl->CastForCell (cell);
            return m_lin2n2ddl;
        case VTK_CELL_TYPE::VTK_QUADRATIC_EDGE:
            m_lin3n3ddl->CastForCell (cell);
            return m_lin3n3ddl;
        case VTK_CELL_TYPE::VTK_TRIANGLE:
            m_tri3n3ddl->CastForCell (cell);
            return m_tri3n3ddl;
        case VTK_CELL_TYPE::VTK_QUADRATIC_TRIANGLE:
            m_tri6n6ddl->CastForCell (cell);
            return m_tri6n6ddl;
        case VTK_CELL_TYPE::VTK_QUAD:
            m_quad4n4ddl->CastForCell (cell);
            return m_quad4n4ddl;
        case VTK_CELL_TYPE::VTK_QUADRATIC_QUAD:
            m_quad8n8ddl->CastForCell (cell);
            return m_quad8n8ddl;
        case VTK_CELL_TYPE::VTK_BIQUADRATIC_QUAD:
            m_quad9n9ddl->CastForCell (cell);
            return m_quad9n9ddl;
        case VTK_CELL_TYPE::VTK_TETRA:
            m_tet4n4ddl->CastForCell (cell);
            return m_tet4n4ddl;
        default:
            m_emp0n0ddl->CastForCell (cell);
            return m_emp0n0ddl;
        }

    }
    else
        return nullptr;
}
