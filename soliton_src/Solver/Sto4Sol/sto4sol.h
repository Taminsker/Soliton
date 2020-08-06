#ifndef STO4SOL_H
#define STO4SOL_H

#include <vector>

class Mesh;

class Sto4Sol
{
public:
    Sto4Sol ();
    ~Sto4Sol ();

    Mesh* mesh;
    std::vector <Mesh *> listobjects;
};

#endif // STO4SOL_H
