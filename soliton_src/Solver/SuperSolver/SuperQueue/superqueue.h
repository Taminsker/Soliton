#ifndef SUPERQUEUE_H
#define SUPERQUEUE_H

#include <vector>
#include <deque>
#include <Algorithms/Math/math.h>
#include <Core/core.h>

class SuperQueue
{
public:
    SuperQueue (Mesh **mesh);
    ~SuperQueue ();
    PlainVector* GetNumber(int num);
    void New (PlainVector& vec);

protected:
    Mesh** m_mesh;
    std::deque<PlainVector*> m_list;
    void CutTheEnd ();
};

#endif // SUPERQUEUE_H
