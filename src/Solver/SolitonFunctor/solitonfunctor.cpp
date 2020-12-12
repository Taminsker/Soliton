#include "solitonfunctor.hpp"

SolitonFunctor NULL_SOLITON_FUNCTOR;

SolitonFunctor::SolitonFunctor ()
{
}

SolitonFunctor::SolitonFunctor (const SolitonFunctor & tocopy) : m_queues (tocopy.m_queues),
                                                                 m_vecs (tocopy.m_vecs)
{
}

SolitonFunctor::~SolitonFunctor ()
{
}

void
SolitonFunctor::LinkTo (std::vector<SolitonQueue *> queues, Mesh * mesh)
{
    m_mesh   = mesh;
    m_queues = queues;
    return;
}

void
SolitonFunctor::LinkTo (std::vector<DenseVector *> vecs, Mesh * mesh)
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
    return Point ();
}

real_t
SolitonFunctor::GetApproximationReal (Point * p, DenseVector * vec, Cell * cell)
{
    if (cell == nullptr || vec == nullptr)
        return 0x0;

    int numPoints = cell->GetNumberOfPoints ();

    // alpha
    DenseVector alpha, values;
    alpha.setZero (numPoints);
    values.setZero (numPoints);

    for (int i = 0; i < numPoints; ++i)
    {
        Point * p_i        = cell->GetPoints ()->at (static_cast<ul_t> (i));
        alpha.coeffRef (i) = EuclidianDist (*p, *p_i);

        values.coeffRef (i) = vec->coeffRef (p_i->GetGlobalIndex ());
    }

    alpha = DenseVector::Ones (numPoints) - alpha / alpha.sum ();

    return (values.cwiseProduct (alpha)).sum ();
}

Point
SolitonFunctor::GetApproximation3DPoint (Point *, DenseVector *, Cell *)
{
    return Point ();
}
