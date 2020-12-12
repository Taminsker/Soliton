#ifndef SRC_SOLVER_SOLITONQUEUE_SOLITONQUEUE_HPP
#define SRC_SOLVER_SOLITONQUEUE_SOLITONQUEUE_HPP

#include <deque>
#include <vector>

#include "../../Algorithms/Math/math.hpp"
#include "../../Core/core.hpp"

class SolitonQueue
{
public:
    SolitonQueue (Mesh ** mesh);
    ~SolitonQueue ();
    DenseVector * GetNumber (int num);
    void          New (DenseVector & vec);

protected:
    Mesh **                   m_mesh;
    std::deque<DenseVector *> m_list;
    void                      CutTheEnd ();
};

#endif /* SRC_SOLVER_SOLITONQUEUE_SOLITONQUEUE_HPP */
