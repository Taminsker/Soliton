#ifndef SUPERQUEUE_H
#define SUPERQUEUE_H

#include <vector>
#include <deque>
#include <Algorithms/Math/math.h>
#include <Core/core.h>

class SolitonQueue
{
public:
    SolitonQueue (Mesh **mesh);
    ~SolitonQueue ();
    PlainVector_eig *GetNumber(int num);
    void New (PlainVector_eig& vec);

protected:
    Mesh** m_mesh;
    std::deque<PlainVector_eig*> m_list;
    void CutTheEnd ();
};

#endif // SUPERQUEUE_H
