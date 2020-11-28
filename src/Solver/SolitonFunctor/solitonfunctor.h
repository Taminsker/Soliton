#ifndef SOLITONFUNCTOR_H
#define SOLITONFUNCTOR_H

#include <Core/core.h>
#include <ProgDef/proddef.h>
#include "../SolitonQueue/solitonqueue.h"

class SolitonFunctor
{
public:
    SolitonFunctor ();
    SolitonFunctor (const SolitonFunctor& tocopy);
    virtual ~SolitonFunctor ();

    void LinkTo (std::vector<SolitonQueue*> queues, Mesh *mesh);
    void LinkTo (std::vector<PlainVector_eig*> vecs, Mesh *mesh);

    virtual real_t ToReal (Point *p, real_t t, Cell *cell);
    virtual Point To3DPoint (Point *p, real_t t, Cell *cell);

protected:
    Mesh                                *m_mesh;
    std::vector<SolitonQueue*>          m_queues;
    std::vector<PlainVector_eig*>       m_vecs;

    real_t GetApproximationReal (Point* p, PlainVector_eig* vec, Cell *cell);
    Point GetApproximation3DPoint (Point* p, PlainVector_eig* vec, Cell *cell);
};

#define SOLITON_FUNCTOR_DEF(X) \
    X () : SolitonFunctor () {} \
    X (const X& tocopy) : SolitonFunctor (tocopy) {} \
    ~X () override {}

#define NULL_SOLITON_FUNCTOR null_soliton_functor
extern SolitonFunctor NULL_SOLITON_FUNCTOR;

#endif // SOLITONFUNCTOR_H
