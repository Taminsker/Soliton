#ifndef DFSTRUCT_H
#define DFSTRUCT_H

#include <vector>
#include <algorithm>

typedef enum
{
    ORDER_1_FORWARD,
    ORDER_2_FORWARD,
    ORDER_3_FORWARD,
    ORDER_4_FORWARD,
    ORDER_5_FORWARD,
    ORDER_6_FORWARD,
    ORDER_7_FORWARD,
    ORDER_8_FORWARD,
    ORDER_1_CENTRAL,
    ORDER_2_CENTRAL,
    ORDER_3_CENTRAL,
    ORDER_4_CENTRAL,
    ORDER_5_CENTRAL,
    ORDER_6_CENTRAL,
    ORDER_7_CENTRAL,
    ORDER_8_CENTRAL,
    ORDER_1_BACKWARD,
    ORDER_2_BACKWARD,
    ORDER_3_BACKWARD,
    ORDER_4_BACKWARD,
    ORDER_5_BACKWARD,
    ORDER_6_BACKWARD,
    ORDER_7_BACKWARD,
    ORDER_8_BACKWARD,
} DF_ORDER_NAMES;

class DFStruct
{

public:
    struct DFIdxCoeff;
    struct DFOrders;
    struct DerivativeStruct;

    DFStruct();
    ~DFStruct ();

    void Get(int der, DF_ORDER_NAMES order, std::vector<int>* idxs, std::vector<double>* coeffs);

private:
    DFStruct (const DFStruct&) = delete;
    DFStruct& operator=(const DFStruct&) = delete;

    DerivativeStruct* Derivative_1;
    DerivativeStruct* Derivative_2;

public:
    struct DFIdxCoeff {
        DFIdxCoeff ();
        ~DFIdxCoeff ();
        DFIdxCoeff (const DFIdxCoeff& tocopy);
        DFIdxCoeff& operator=(const DFIdxCoeff& tocopy);

        std::vector<int>* idxs;
        std::vector<double>* coeffs;
    };

    struct DFOrders {
        DFOrders ();
        ~DFOrders ();

        DFIdxCoeff* Order1;
        DFIdxCoeff* Order2;
        DFIdxCoeff* Order3;
        DFIdxCoeff* Order4;
        DFIdxCoeff* Order5;
        DFIdxCoeff* Order6;
        DFIdxCoeff* Order7;
        DFIdxCoeff* Order8;
    };

    struct DerivativeStruct {
        DerivativeStruct ();
        ~DerivativeStruct ();

        DFOrders* Backward;
        DFOrders* Central;
        DFOrders* Forward;
    };
};

using DFPair = DFStruct::DFIdxCoeff;

void DFOrderBuild (int der, DF_ORDER_NAMES order, std::vector<int>* idxs, std::vector<double>* coeffs);

void Reverse (DFStruct::DFOrders* O, DFStruct::DFOrders* R);

void Reverse (DFStruct::DFIdxCoeff* D, DFStruct::DFIdxCoeff* R);

template <typename T>
std::vector<T> operator* (T value, std::vector<T>& vec)
{
    std::vector<T> R(vec.size (), T(0));

    for (std::size_t i = 0; i < vec.size (); ++i)
    {
        R.at (i) = value * vec.at (i);
    }

    return R;
}


#endif // DFSTRUCT_H
