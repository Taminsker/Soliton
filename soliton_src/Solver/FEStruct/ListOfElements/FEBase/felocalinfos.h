#ifndef FELOCALINFOS_H
#define FELOCALINFOS_H

#include <vector>
#include <Algorithms/Math/math.h>
#include <Core/core.h>

class FELocalInfos
{
public:
    FELocalInfos ();
    ~FELocalInfos ();

    Matrix3x3 Jac;
    Matrix3x3 JacInvT;
    Matrix3x3 Jac_cmp;
    Matrix3x3 JacInvT_cmp;
    double detJac;
    Point normalEdge;
    Point tangentEdge;
    Point normalCell;
    Point normalEdge_ref;
    Point tangentEdge_ref;
    Point normalCell_ref;
};

#endif // FELOCALINFOS_H
