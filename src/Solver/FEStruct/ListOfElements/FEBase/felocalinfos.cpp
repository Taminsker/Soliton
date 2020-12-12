#include "felocalinfos.hpp"

FELocalInfos::FELocalInfos () : Jac (Matrix3x3::Identity ()),
                                JacInvT (Matrix3x3::Identity ()),
                                Jac_cmp (Matrix3x3::Identity ()),
                                JacInvT_cmp (Matrix3x3::Identity ()),
                                detJac (0.0),
                                normalEdge ({}),
                                tangentEdge ({}),
                                normalCell ({}),
                                normalEdge_ref ({}),
                                tangentEdge_ref ({}),
                                normalCell_ref ({})
{
}

FELocalInfos::~FELocalInfos ()
{
}
