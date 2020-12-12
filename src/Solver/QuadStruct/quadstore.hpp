#ifndef SRC_SOLVER_QUADSTRUCT_QUADSTORE_HPP
#define SRC_SOLVER_QUADSTRUCT_QUADSTORE_HPP

#include <vector>

#include "../../Core/core.hpp"

class FEBase;

class QuadStore
{
public:
    struct QuadObject;

    QuadStore ();
    ~QuadStore ();
    QuadObject * Get (FEBase * object);

private:
    QuadStore (const QuadStore &) = delete;
    QuadStore & operator= (const QuadStore &) = delete;

    QuadObject *              m_quad_empty;
    QuadObject *              m_quad_vertex;
    std::vector<QuadObject *> m_quad_line;
    std::vector<QuadObject *> m_quad_triangle;
    std::vector<QuadObject *> m_quad_quadrangle;

public:
    struct QuadObject
    {
        QuadObject (PHYSICAL_CELL_TYPE _type = PHYSICAL_CELL_TYPE::EMPTY, ul_t _npts = 0, ul_t _order = 0);
        ~QuadObject ();

        void Print ();

        PHYSICAL_CELL_TYPE  type;
        ul_t                npts;
        ul_t                order;
        std::vector<real_t> w;
        std::vector<Point>  pts;
    };
};

QuadStore::QuadObject * QuadrangleFromLine (QuadStore::QuadObject * line);
#endif /* SRC_SOLVER_QUADSTRUCT_QUADSTORE_HPP */
