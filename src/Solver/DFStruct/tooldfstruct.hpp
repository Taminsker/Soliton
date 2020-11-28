#ifndef TOOLDFSTRUCT_H
#define TOOLDFSTRUCT_H

#include <vector>
#include <algorithm>

#include <Solver/DFStruct/dfstore.h>

void DFOrderBuild (int der, DF_ORDER_NAMES order, std::vector<int>* idxs, std::vector<real_t>* coeffs);

void Reverse (DFStore::DFOrders& O, DFStore::DFOrders& R);

void Reverse (DFStore::DFIdxCoeff& D, DFStore::DFIdxCoeff& R);

template <typename T>
std::vector<T> operator* (T value, const std::vector<T>& vec)
{
    std::vector<T> R(vec.size (), T(0));

    for (ul_t i = 0; i < vec.size (); ++i)
    {
        R.at (i) = value * vec.at (i);
    }

    return R;
}

#endif // TOOLDFSTRUCT_H
