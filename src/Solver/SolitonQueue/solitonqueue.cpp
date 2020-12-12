#include "solitonqueue.hpp"

SolitonQueue::SolitonQueue (Mesh ** mesh) : m_mesh (mesh),
                                            m_list ({})
{
    int numPoints = (*m_mesh)->GetNumberOfPoints ();

    for (ul_t id = 0; id < SOLITONQUEUE_MAX; id++)
    {
        DenseVector * obj = new DenseVector ();
        obj->resize (numPoints);
        obj->setZero ();
        m_list.push_front (obj);
    }
}

SolitonQueue::~SolitonQueue ()
{
    for (DenseVector * obj : m_list)
        delete obj;

    m_list.clear ();
}

DenseVector *
SolitonQueue::GetNumber (int num)
{
    return m_list.at (ul_t (std::abs (num)));
}

void
SolitonQueue::New (DenseVector & vec)
{
    m_list.push_front (new DenseVector (vec));
    return CutTheEnd ();
}

void
SolitonQueue::CutTheEnd ()
{
    while (m_list.size () != SOLITONQUEUE_MAX)
    {
        delete m_list.back ();
        m_list.erase (m_list.end ());
    }

    return;
}
