#include <Solver/SolitonFE/SolitonFEItem/solitonfeitem.h>
#include <Solver/SolitonFE/SolitonFEContainer/solitonfecontainer.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>
#include <Algorithms/algorithms.h>

template<>
void
SolitonFEItem<ITEM_T::EMPTY>::Compute (real_t)
{
    return;
}
