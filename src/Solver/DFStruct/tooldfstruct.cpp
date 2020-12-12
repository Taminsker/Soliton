#include "tooldfstruct.hpp"

void
Reverse (DFStore::DFOrders & O, DFStore::DFOrders & R)
{
    Reverse (O.Order1, R.Order1);
    Reverse (O.Order2, R.Order2);
    Reverse (O.Order3, R.Order3);
    Reverse (O.Order4, R.Order4);
    Reverse (O.Order5, R.Order5);
    Reverse (O.Order6, R.Order6);
    Reverse (O.Order7, R.Order7);
    Reverse (O.Order8, R.Order8);

    return;
}

void
Reverse (DFStore::DFIdxCoeff & D, DFStore::DFIdxCoeff & R)
{
    R = D;

    R.idxs   = -1 * R.idxs;
    R.coeffs = -1. * R.coeffs;

    std::reverse (R.idxs.begin (), R.idxs.end ());
    std::reverse (R.coeffs.begin (), R.coeffs.end ());

    return;
}

void
DFStore::Get (int der, DF_ORDER_NAMES order, std::vector<int> * idxs, std::vector<real_t> * coeffs)
{
    if (der == 1)
    {
        switch (order)
        {
            case ORDER_1_FORWARD:
                *idxs   = Derivative_1.Forward.Order1.idxs;
                *coeffs = Derivative_1.Forward.Order1.coeffs;
                break;
            case ORDER_1_CENTRAL:
                *idxs   = Derivative_1.Central.Order1.idxs;
                *coeffs = Derivative_1.Central.Order1.coeffs;
                break;
            case ORDER_1_BACKWARD:
                *idxs   = Derivative_1.Backward.Order1.idxs;
                *coeffs = Derivative_1.Backward.Order1.coeffs;
                break;
            case ORDER_2_FORWARD:
                *idxs   = Derivative_1.Forward.Order2.idxs;
                *coeffs = Derivative_1.Forward.Order2.coeffs;
                break;
            case ORDER_2_CENTRAL:
                *idxs   = Derivative_1.Central.Order2.idxs;
                *coeffs = Derivative_1.Central.Order2.coeffs;
                break;
            case ORDER_2_BACKWARD:
                *idxs   = Derivative_1.Backward.Order2.idxs;
                *coeffs = Derivative_1.Backward.Order2.coeffs;
                break;
            case ORDER_3_FORWARD:
                *idxs   = Derivative_1.Forward.Order3.idxs;
                *coeffs = Derivative_1.Forward.Order3.coeffs;
                break;
            case ORDER_3_CENTRAL:
                *idxs   = Derivative_1.Central.Order3.idxs;
                *coeffs = Derivative_1.Central.Order3.coeffs;
                break;
            case ORDER_3_BACKWARD:
                *idxs   = Derivative_1.Backward.Order3.idxs;
                *coeffs = Derivative_1.Backward.Order3.coeffs;
                break;
            case ORDER_4_FORWARD:
                *idxs   = Derivative_1.Forward.Order4.idxs;
                *coeffs = Derivative_1.Forward.Order4.coeffs;
                break;
            case ORDER_4_CENTRAL:
                *idxs   = Derivative_1.Central.Order4.idxs;
                *coeffs = Derivative_1.Central.Order4.coeffs;
                break;
            case ORDER_4_BACKWARD:
                *idxs   = Derivative_1.Backward.Order4.idxs;
                *coeffs = Derivative_1.Backward.Order4.coeffs;
                break;
            case ORDER_5_FORWARD:
                *idxs   = Derivative_1.Forward.Order5.idxs;
                *coeffs = Derivative_1.Forward.Order5.coeffs;
                break;
            case ORDER_5_CENTRAL:
                *idxs   = Derivative_1.Central.Order5.idxs;
                *coeffs = Derivative_1.Central.Order5.coeffs;
                break;
            case ORDER_5_BACKWARD:
                *idxs   = Derivative_1.Backward.Order5.idxs;
                *coeffs = Derivative_1.Backward.Order5.coeffs;
                break;
            case ORDER_6_FORWARD:
                *idxs   = Derivative_1.Forward.Order6.idxs;
                *coeffs = Derivative_1.Forward.Order6.coeffs;
                break;
            case ORDER_6_CENTRAL:
                *idxs   = Derivative_1.Central.Order6.idxs;
                *coeffs = Derivative_1.Central.Order6.coeffs;
                break;
            case ORDER_6_BACKWARD:
                *idxs   = Derivative_1.Backward.Order6.idxs;
                *coeffs = Derivative_1.Backward.Order6.coeffs;
                break;
            case ORDER_7_FORWARD:
                *idxs   = Derivative_1.Forward.Order7.idxs;
                *coeffs = Derivative_1.Forward.Order7.coeffs;
                break;
            case ORDER_7_CENTRAL:
                *idxs   = Derivative_1.Central.Order7.idxs;
                *coeffs = Derivative_1.Central.Order7.coeffs;
                break;
            case ORDER_7_BACKWARD:
                *idxs   = Derivative_1.Backward.Order7.idxs;
                *coeffs = Derivative_1.Backward.Order7.coeffs;
                break;
            case ORDER_8_FORWARD:
                *idxs   = Derivative_1.Forward.Order8.idxs;
                *coeffs = Derivative_1.Forward.Order8.coeffs;
                break;
            case ORDER_8_CENTRAL:
                *idxs   = Derivative_1.Central.Order8.idxs;
                *coeffs = Derivative_1.Central.Order8.coeffs;
                break;
            case ORDER_8_BACKWARD:
                *idxs   = Derivative_1.Backward.Order8.idxs;
                *coeffs = Derivative_1.Backward.Order8.coeffs;
                break;
        }
    }
    else if (der == 2)
    {
        switch (order)
        {
            case ORDER_1_FORWARD:
                *idxs   = Derivative_2.Forward.Order1.idxs;
                *coeffs = Derivative_2.Forward.Order1.coeffs;
                break;
            case ORDER_1_CENTRAL:
                *idxs   = Derivative_2.Central.Order1.idxs;
                *coeffs = Derivative_2.Central.Order1.coeffs;
                break;
            case ORDER_1_BACKWARD:
                *idxs   = Derivative_2.Backward.Order1.idxs;
                *coeffs = Derivative_2.Backward.Order1.coeffs;
                break;
            case ORDER_2_FORWARD:
                *idxs   = Derivative_2.Forward.Order2.idxs;
                *coeffs = Derivative_2.Forward.Order2.coeffs;
                break;
            case ORDER_2_CENTRAL:
                *idxs   = Derivative_2.Central.Order2.idxs;
                *coeffs = Derivative_2.Central.Order2.coeffs;
                break;
            case ORDER_2_BACKWARD:
                *idxs   = Derivative_2.Backward.Order2.idxs;
                *coeffs = Derivative_2.Backward.Order2.coeffs;
                break;
            case ORDER_3_FORWARD:
                *idxs   = Derivative_2.Forward.Order3.idxs;
                *coeffs = Derivative_2.Forward.Order3.coeffs;
                break;
            case ORDER_3_CENTRAL:
                *idxs   = Derivative_2.Central.Order3.idxs;
                *coeffs = Derivative_2.Central.Order3.coeffs;
                break;
            case ORDER_3_BACKWARD:
                *idxs   = Derivative_2.Backward.Order3.idxs;
                *coeffs = Derivative_2.Backward.Order3.coeffs;
                break;
            case ORDER_4_FORWARD:
                *idxs   = Derivative_2.Forward.Order4.idxs;
                *coeffs = Derivative_2.Forward.Order4.coeffs;
                break;
            case ORDER_4_CENTRAL:
                *idxs   = Derivative_2.Central.Order4.idxs;
                *coeffs = Derivative_2.Central.Order4.coeffs;
                break;
            case ORDER_4_BACKWARD:
                *idxs   = Derivative_2.Backward.Order4.idxs;
                *coeffs = Derivative_2.Backward.Order4.coeffs;
                break;
            case ORDER_5_FORWARD:
                *idxs   = Derivative_2.Forward.Order5.idxs;
                *coeffs = Derivative_2.Forward.Order5.coeffs;
                break;
            case ORDER_5_CENTRAL:
                *idxs   = Derivative_2.Central.Order5.idxs;
                *coeffs = Derivative_2.Central.Order5.coeffs;
                break;
            case ORDER_5_BACKWARD:
                *idxs   = Derivative_2.Backward.Order5.idxs;
                *coeffs = Derivative_2.Backward.Order5.coeffs;
                break;
            case ORDER_6_FORWARD:
                *idxs   = Derivative_2.Forward.Order6.idxs;
                *coeffs = Derivative_2.Forward.Order6.coeffs;
                break;
            case ORDER_6_CENTRAL:
                *idxs   = Derivative_2.Central.Order6.idxs;
                *coeffs = Derivative_2.Central.Order6.coeffs;
                break;
            case ORDER_6_BACKWARD:
                *idxs   = Derivative_2.Backward.Order6.idxs;
                *coeffs = Derivative_2.Backward.Order6.coeffs;
                break;
            case ORDER_7_FORWARD:
                *idxs   = Derivative_2.Forward.Order7.idxs;
                *coeffs = Derivative_2.Forward.Order7.coeffs;
                break;
            case ORDER_7_CENTRAL:
                *idxs   = Derivative_2.Central.Order7.idxs;
                *coeffs = Derivative_2.Central.Order7.coeffs;
                break;
            case ORDER_7_BACKWARD:
                *idxs   = Derivative_2.Backward.Order7.idxs;
                *coeffs = Derivative_2.Backward.Order7.coeffs;
                break;
            case ORDER_8_FORWARD:
                *idxs   = Derivative_2.Forward.Order8.idxs;
                *coeffs = Derivative_2.Forward.Order8.coeffs;
                break;
            case ORDER_8_CENTRAL:
                *idxs   = Derivative_2.Central.Order8.idxs;
                *coeffs = Derivative_2.Central.Order8.coeffs;
                break;
            case ORDER_8_BACKWARD:
                *idxs   = Derivative_2.Backward.Order8.idxs;
                *coeffs = Derivative_2.Backward.Order8.coeffs;
                break;
        }
    }

    return;
}
