#include <Solver/SuperSolver/SuperItem/superitem.h>
#include <Solver/FEStruct/festruct.h>
#include <Solver/QuadStruct/quadstruct.h>
#include <Algorithms/algorithms.h>


template<>
void
SuperItem::InternalCompute<ITEM_T::EMPTY> ()
{
    SuperItem *itemThis = this;
    VOID_USE(itemThis);

    return;
}
