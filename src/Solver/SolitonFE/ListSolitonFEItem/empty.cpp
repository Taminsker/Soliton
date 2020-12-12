#include "../../../Algorithms/algorithms.hpp"
#include "../../FEStruct/festruct.hpp"
#include "../../QuadStruct/quadstruct.hpp"
#include "../SolitonFEContainer/solitonfecontainer.hpp"
#include "../SolitonFEItem/solitonfeitem.hpp"

template <>
void
    SolitonFEItem<ITEM_T::EMPTY>::Compute (real_t)
{
    return;
}
