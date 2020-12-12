#ifndef SRC_SOLVER_MESHSTORAGE_MESHSTORAGE_HPP
#define SRC_SOLVER_MESHSTORAGE_MESHSTORAGE_HPP

#include <vector>

#include "../../solitonheader.hpp"

class Mesh;
class InputDatStruct;

class MeshStorage
{
public:
    MeshStorage ();
    ~MeshStorage ();

    void                  PushBack (Mesh * object);
    std::vector<Mesh *> * GetListOfObjects ();
    Mesh *                GetMainMesh ();
    Mesh *                GetMeshObjectAt (ul_t id);
    Mesh **               GetMainMesh_inc ();
    Mesh **               GetObjectAt_inc (ul_t id);

protected:
    Mesh *              m_principal_mesh;
    std::vector<Mesh *> m_list_objects;
};

#endif /* SRC_SOLVER_MESHSTORAGE_MESHSTORAGE_HPP */
