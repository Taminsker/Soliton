#include "sto4sol.h"
#include <Core/core.h>

Sto4Sol::Sto4Sol() :
    mesh (new Mesh()),
    listobjects ({})
{
    mesh->SetName ("mesh");
}

Sto4Sol::~Sto4Sol ()
{
    delete mesh;
    for (Mesh* obj : listobjects)
        delete obj;
    listobjects.clear ();
}
