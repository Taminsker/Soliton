#ifndef FESTORE_H
#define FESTORE_H

#include <Core/core.h>
#include <Solver/FEStruct/ListOfElements/listofelements.h>

class FEStore
{
public:
    FEStore ();
    ~FEStore ();

    FEBase* GetElementFor (Cell* cell, FE_CLASS_TYPE typeOfEle = FE_CLASS_TYPE::LAGRANGE);

private:
    Emp0N0DDL*   m_emp0n0ddl;
    Ver1N1DDL*   m_ver1n1ddl;
    Lin2N2DDL*   m_lin2n2ddl;
    Lin3N3DDL*   m_lin3n3ddl;
    Tri3N3DDL*   m_tri3n3ddl;
    Tri6N6DDL*   m_tri6n6ddl;
    Quad4N4DDL*   m_quad4n4ddl;
    Quad8N8DDL*   m_quad8n8ddl;
    Quad9N9DDL*   m_quad9n9ddl;
    Tet4N4DDL*   m_tet4n4ddl;
};

#endif // FESTORE_H
