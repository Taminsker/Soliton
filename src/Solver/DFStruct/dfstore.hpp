#ifndef SRC_SOLVER_DFSTRUCT_DFSTORE_HPP
#define SRC_SOLVER_DFSTRUCT_DFSTORE_HPP

#include <algorithm>
#include <vector>

#include "../../solitonheader.hpp"

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

class DFStore
{
public:
    class DFIdxCoeff
    {
    public:
        DFIdxCoeff ();
        ~DFIdxCoeff ();
        DFIdxCoeff (const DFIdxCoeff & tocopy);
        DFIdxCoeff & operator= (const DFIdxCoeff & tocopy);

        std::vector<int>    idxs;
        std::vector<real_t> coeffs;
    };

    class DFOrders
    {
    public:
        DFOrders ();
        ~DFOrders ();

        DFIdxCoeff Order1;
        DFIdxCoeff Order2;
        DFIdxCoeff Order3;
        DFIdxCoeff Order4;
        DFIdxCoeff Order5;
        DFIdxCoeff Order6;
        DFIdxCoeff Order7;
        DFIdxCoeff Order8;
    };

    class DerivativeStruct
    {
    public:
        DerivativeStruct ();
        ~DerivativeStruct ();

        DFOrders Backward;
        DFOrders Central;
        DFOrders Forward;
    };

    DFStore ();
    ~DFStore ();

    void Get (int der, DF_ORDER_NAMES order, std::vector<int> * idxs, std::vector<real_t> * coeffs);

private:
    DFStore (const DFStore &) = delete;
    DFStore & operator= (const DFStore &) = delete;

    DerivativeStruct Derivative_1;
    DerivativeStruct Derivative_2;
};

#endif /* SRC_SOLVER_DFSTRUCT_DFSTORE_HPP */
