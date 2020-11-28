#include "solitonfunctor.h"

SolitonFunctor NULL_SOLITON_FUNCTOR;

SolitonFunctor::SolitonFunctor()
{}

SolitonFunctor::SolitonFunctor (const SolitonFunctor &tocopy) :
    m_queues (tocopy.m_queues),
    m_vecs (tocopy.m_vecs)
{}

SolitonFunctor::~SolitonFunctor ()
{}

void SolitonFunctor::LinkTo(std::vector<SolitonQueue *> queues, Mesh *mesh)
{
    m_mesh = mesh;
    m_queues = queues;
    return;
}

void SolitonFunctor::LinkTo(std::vector<PlainVector_eig *> vecs, Mesh *mesh)
{
    m_mesh = mesh;
    m_vecs = vecs;
    return;
}

real_t
SolitonFunctor::ToReal (Point *, real_t, Cell *)
{
    return 0x0;
}

Point
SolitonFunctor::To3DPoint (Point *, real_t, Cell *)
{
    return Point();
}

real_t
SolitonFunctor::GetApproximationReal (Point* p, PlainVector_eig *vec, Cell *cell)
{
    if (cell == nullptr || vec == nullptr)
        return 0x0;

    int numPoints = cell->GetNumberOfPoints ();

    // alpha
    PlainVector_eig alpha, values;
    alpha.setZero (numPoints);
    values.setZero (numPoints);

    for (int i = 0; i < numPoints; ++i)
    {
        Point *p_i = cell->GetPoints ()->at (static_cast<ul_t>(i));
        alpha.coeffRef (i) =  EuclidianDist (*p, *p_i);

        values.coeffRef (i) = vec->coeffRef (p_i->GetGlobalIndex ());
    }

    alpha = PlainVector_eig::Ones (numPoints) - alpha / alpha.sum ();

    return (values.cwiseProduct (alpha)).sum ();
}

Point
SolitonFunctor::GetApproximation3DPoint (Point *, PlainVector_eig *, Cell *)
{
    return Point ();
}
