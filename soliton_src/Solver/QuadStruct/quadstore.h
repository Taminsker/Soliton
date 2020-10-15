#ifndef QUADSTORE_H
#define QUADSTORE_H

#include <vector>

#include <Core/core.h>
#include <Solver/FEStruct/ListOfElements/listofelements.h>

class FEBase;

class QuadStore
{
public:
   struct QuadObject;

    QuadStore ();
    ~QuadStore ();
    QuadObject* Get (FEBase* object);

private:
    QuadStore (const QuadStore&) = delete;
    QuadStore& operator=(const QuadStore&) = delete;

    QuadObject* m_quad_empty;
    QuadObject* m_quad_vertex;
    std::vector<QuadObject*> m_quad_line;
    std::vector<QuadObject*> m_quad_triangle;
    std::vector<QuadObject*> m_quad_quadrangle;

public:
    struct QuadObject
    {
      QuadObject (PHYSICAL_CELL_TYPE _type = PHYSICAL_CELL_TYPE::EMPTY, std::size_t _npts = 0, std::size_t _order = 0);
      ~QuadObject ();

      void Print();

      PHYSICAL_CELL_TYPE type;
      std::size_t npts;
      std::size_t order;
      std::vector <double> w;
      std::vector<Point> pts;
    };
};

QuadStore::QuadObject* QuadrangleFromLine(QuadStore::QuadObject* line);
#endif // QUADSTORE_H
