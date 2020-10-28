#include "dfstore.h"
#include <Solver/DFStruct/tooldfstruct.h>

DFStore::DFIdxCoeff::DFIdxCoeff () :
    idxs ({}),
    coeffs ({})
{}

DFStore::DFIdxCoeff::~DFIdxCoeff ()
{}

DFStore::DFIdxCoeff::DFIdxCoeff (const DFStore::DFIdxCoeff& tocopy) :
    idxs (tocopy.idxs),
    coeffs (tocopy.coeffs)
{}

DFStore::DFIdxCoeff& DFStore::DFIdxCoeff::operator=(const DFStore::DFIdxCoeff& tocopy)
{
    idxs = tocopy.idxs;
    coeffs = tocopy.coeffs;

    return *this;
}


DFStore::DFOrders::DFOrders () :
    Order1 (),
    Order2 (),
    Order3 (),
    Order4 (),
    Order5 (),
    Order6 (),
    Order7 (),
    Order8 ()
{}

DFStore::DFOrders::~DFOrders ()
{}

DFStore::DerivativeStruct::DerivativeStruct () :
    Backward (),
    Central  (),
    Forward  ()
{}

DFStore::DerivativeStruct::~DerivativeStruct ()
{}


DFStore::DFStore () :
    Derivative_1 (),
    Derivative_2 ()
{
    // Derivée première
    // Forward
    Derivative_1.Forward.Order1.idxs = {0, 1};
    Derivative_1.Forward.Order1.coeffs = {-1., 1.};

    Derivative_1.Forward.Order2.idxs = {0, 1, 2};
    Derivative_1.Forward.Order2.coeffs = {-3./2., 2., -1./2.};

    Derivative_1.Forward.Order3.idxs = {0, 1, 2, 3};
    Derivative_1.Forward.Order3.coeffs = {-11./6., 3., -3./2., 1./3.};

    Derivative_1.Forward.Order4.idxs = {0, 1, 2, 3, 4};
    Derivative_1.Forward.Order4.coeffs = {-25./12., 4., -3., 4./3., -1./4.};

    Derivative_1.Forward.Order5.idxs = {0, 1, 2, 3, 4, 5};
    Derivative_1.Forward.Order5.coeffs = {-137./60., 5., -5., 10./3., -5./4., 1./5.};

    Derivative_1.Forward.Order6.idxs = {0, 1, 2, 3, 4, 5, 6};
    Derivative_1.Forward.Order6.coeffs = {-49./20., 6., -15./2., 20./3., -15./4., 6./5., -1./6.};

    // Central
    Derivative_1.Central.Order2.idxs = {-1, 1};
    Derivative_1.Central.Order2.coeffs = {-1./2., 1./2.};

    Derivative_1.Central.Order4.idxs = {-2, -1, 1, 2};
    Derivative_1.Central.Order4.coeffs = {1./12., -2./3., 2./3., -1./12.};

    Derivative_1.Central.Order6.idxs = {-3, -2, -1, 1, 2, 3};
    Derivative_1.Central.Order6.coeffs = {-1./60., 3./20., -3./4., 3./4., -3./20., 1./60.};

    Derivative_1.Central.Order8.idxs = {-4, -3, -2, -1, 1, 2, 3, 4};
    Derivative_1.Central.Order8.coeffs = {1./280., -4./105., 1./5., -4./5., 4./5., -1./5., 4./105., -1./280.};

    // Backward
    Reverse (Derivative_1.Forward, Derivative_1.Backward);


    // Derivée seconde
    // Forward
    Derivative_2.Forward.Order1.idxs = {0, 1, 2};
    Derivative_2.Forward.Order1.coeffs = {1., 2., -1./2.};

    Derivative_2.Forward.Order2.idxs = {0, 1, 2, 3};
    Derivative_2.Forward.Order2.coeffs = {2., -5., 4., -1.};

    Derivative_2.Forward.Order3.idxs = {0, 1, 2, 3, 4};
    Derivative_2.Forward.Order3.coeffs = {35./12., -26./3., 19./2., -14./3., 11./12.};

    Derivative_2.Forward.Order4.idxs = {0, 1, 2, 3, 4, 5};
    Derivative_2.Forward.Order4.coeffs = {15./4., -77./6., 107./6., -13., 61./12., -5./6.};

    Derivative_2.Forward.Order5.idxs = {0, 1, 2, 3, 4, 5, 6};
    Derivative_2.Forward.Order5.coeffs = {203./45., -87./5., 117./4., -254./9., 33./2., -27./5., 137./180.};

    Derivative_2.Forward.Order6.idxs = {0, 1, 2, 3, 4, 5, 6, 7};
    Derivative_2.Forward.Order6.coeffs = {469./90., -223./10., 879./20., -949./18., 41., -201./10., 1019./180., -7./10.};

    // Central
    Derivative_2.Central.Order2.idxs = {-1, 0, 1};
    Derivative_2.Central.Order2.coeffs = {1., -2.,  1.};

    Derivative_2.Central.Order4.idxs = {-2, -1, 0, 1, 2};
    Derivative_2.Central.Order4.coeffs = {-1./12., 4./3., -5./2., 4./3., -1./12.};

    Derivative_2.Central.Order6.idxs = {-3, -2, -1, 0, 1, 2, 3};
    Derivative_2.Central.Order6.coeffs = {1./90., -3./20., 3./2., -49./18., 3./2., -3./20., 1./90.};

    Derivative_2.Central.Order8.idxs = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
    Derivative_2.Central.Order8.coeffs = {-1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560.};

    // Backward
    Reverse (Derivative_2.Forward, Derivative_2.Backward);
}

DFStore::~DFStore ()
{}

