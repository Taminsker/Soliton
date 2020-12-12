#ifndef SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEBASE_FELOCALINFOS_HPP
#define SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEBASE_FELOCALINFOS_HPP

#include <vector>

#include "../../../../Algorithms/Math/math.hpp"
#include "../../../../Core/core.hpp"

class FELocalInfos
{
public:
    FELocalInfos ();
    ~FELocalInfos ();

    Matrix3x3 Jac;
    Matrix3x3 JacInvT;
    Matrix3x3 Jac_cmp;
    Matrix3x3 JacInvT_cmp;
    real_t    detJac;
    Point     normalEdge;
    Point     tangentEdge;
    Point     normalCell;
    Point     normalEdge_ref;
    Point     tangentEdge_ref;
    Point     normalCell_ref;
};

#endif /* SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEBASE_FELOCALINFOS_HPP */
