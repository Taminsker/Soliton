#ifndef MESHSTORAGE_H
#define MESHSTORAGE_H

#include <ProgDef/proddef.h>
#include <vector>

class Mesh;
class InputDatStruct;

class MeshStorage
{
public:
    MeshStorage ();
    ~MeshStorage ();

    void PushBack (Mesh * object);
    std::vector<Mesh *>* GetListOfObjects ();
    Mesh *GetMainMesh ();
    Mesh *GetMeshObjectAt (ul_t id);
    Mesh **GetMainMesh_inc ();
    Mesh **GetObjectAt_inc (ul_t id);

protected:
    Mesh* m_principal_mesh;
    std::vector <Mesh *> m_list_objects;
};

#endif // MESHSTORAGE_H
