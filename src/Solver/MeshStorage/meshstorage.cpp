#include "meshstorage.hpp"

#include "../../Core/core.hpp"

MeshStorage::MeshStorage () : m_principal_mesh (new Mesh ()),
                              m_list_objects ({})
{
    GetMainMesh ()->SetName ("mesh");
}

MeshStorage::~MeshStorage ()
{
    delete m_principal_mesh;

    for (Mesh * obj : m_list_objects)
        delete obj;
    m_list_objects.clear ();
}

void
MeshStorage::PushBack (Mesh * mesh)
{
    m_list_objects.push_back (mesh);
    return;
}

Mesh *
MeshStorage::GetMainMesh ()
{
    return m_principal_mesh;
}

std::vector<Mesh *> *
MeshStorage::GetListOfObjects ()
{
    return &m_list_objects;
}

Mesh *
MeshStorage::GetMeshObjectAt (ul_t id)
{
    return m_list_objects.at (id);
}

Mesh **
MeshStorage::GetMainMesh_inc ()
{
    return &m_principal_mesh;
}
Mesh **
MeshStorage::GetObjectAt_inc (ul_t id)
{
    return &m_list_objects.at (id);
}
