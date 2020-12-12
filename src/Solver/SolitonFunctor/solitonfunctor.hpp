#ifndef SRC_SOLVER_SOLITONFUNCTOR_SOLITONFUNCTOR_HPP
#define SRC_SOLVER_SOLITONFUNCTOR_SOLITONFUNCTOR_HPP

#include "../../Core/core.hpp"
#include "../../solitonheader.hpp"
#include "../SolitonQueue/solitonqueue.hpp"

class SolitonFunctor
{
public:
    SolitonFunctor ();
    SolitonFunctor (const SolitonFunctor & tocopy);
    virtual ~SolitonFunctor ();

    void LinkTo (std::vector<SolitonQueue *> queues, Mesh * mesh);
    void LinkTo (std::vector<DenseVector *> vecs, Mesh * mesh);

    virtual real_t ToReal (Point * p, real_t t, Cell * cell);
    virtual Point  To3DPoint (Point * p, real_t t, Cell * cell);

protected:
    Mesh *                      m_mesh;
    std::vector<SolitonQueue *> m_queues;
    std::vector<DenseVector *>  m_vecs;

    real_t GetApproximationReal (Point * p, DenseVector * vec, Cell * cell);
    Point  GetApproximation3DPoint (Point * p, DenseVector * vec, Cell * cell);
};

#define SOLITON_FUNCTOR_DEF(X)                        \
    X () : SolitonFunctor () {}                       \
    X (const X & tocopy) : SolitonFunctor (tocopy) {} \
    ~X () override {}

#define NULL_SOLITON_FUNCTOR null_soliton_functor
extern SolitonFunctor NULL_SOLITON_FUNCTOR;

#endif /* SRC_SOLVER_SOLITONFUNCTOR_SOLITONFUNCTOR_HPP */
