//#include "queuesolutions.h"

//QueueSolutions::QueueSolutions (Mesh **mesh) :
//    m_mesh (mesh),
//    m_list ({})
//{

//    int numPoints = (*m_mesh)->GetNumberOfPoints ();

//    for (std::size_t id = 0; id < SOLQUEUE_MAX; id++)
//    {
//        PlainVector* obj = new PlainVector ();
//        obj->resize(numPoints);
//        obj->setZero ();
//        m_list.push_front (obj);
//    }
//}

//QueueSolutions::~QueueSolutions ()
//{
//    for (PlainVector* obj : m_list)
//        delete obj;

//    m_list.clear ();
//}

//PlainVector* QueueSolutions::GetNumber (int num)
//{
//    return m_list.at (std::size_t (std::abs(num)));
//}

//void QueueSolutions::New (PlainVector& vec)
//{
//    m_list.push_front (new PlainVector (vec));
//    return CutTheEnd ();
//}

//void QueueSolutions::CutTheEnd ()
//{
//    while (m_list.size () != SOLQUEUE_MAX)
//    {
//        delete m_list.back ();
//        m_list.erase (m_list.end ());
//    }

//    return;
//}
