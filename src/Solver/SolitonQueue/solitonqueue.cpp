#include "solitonqueue.h"

SolitonQueue::SolitonQueue (Mesh **mesh) :
    m_mesh (mesh),
    m_list ({})
{

    int numPoints = (*m_mesh)->GetNumberOfPoints ();

    for (ul_t id = 0; id < SOLITONQUEUE_MAX; id++)
    {
        PlainVector_eig *obj = new PlainVector_eig ();
        obj->resize(numPoints);
        obj->setZero ();
        m_list.push_front (obj);
    }
}

SolitonQueue::~SolitonQueue ()
{
    for (PlainVector_eig *obj : m_list)
        delete obj;

    m_list.clear ();
}

PlainVector_eig *SolitonQueue::GetNumber (int num)
{
    return m_list.at (ul_t (std::abs(num)));
}

void SolitonQueue::New (PlainVector_eig &vec)
{
    m_list.push_front (new PlainVector_eig (vec));
    return CutTheEnd ();
}

void SolitonQueue::CutTheEnd ()
{
    while (m_list.size () != SOLITONQUEUE_MAX)
    {
        delete m_list.back ();
        m_list.erase (m_list.end ());
    }

    return;
}
