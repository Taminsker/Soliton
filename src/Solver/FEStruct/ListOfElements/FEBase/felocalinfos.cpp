#include "felocalinfos.h"

FELocalInfos::FELocalInfos () :
    Jac (Matrix3x3_eig::Identity ()),
    JacInvT (Matrix3x3_eig::Identity ()),
    Jac_cmp (Matrix3x3_eig::Identity ()),
    JacInvT_cmp (Matrix3x3_eig::Identity ()),
    detJac (0.0),
    normalEdge ({}),
    tangentEdge ({}),
    normalCell ({}),
    normalEdge_ref ({}),
    tangentEdge_ref ({}),
    normalCell_ref ({})
{}

FELocalInfos::~FELocalInfos ()
{}


