#ifndef SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEHERMITE_FEHERMITE_HPP
#define SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEHERMITE_FEHERMITE_HPP

#include "../PhysicalElements/physicalelements.hpp"

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
class FEHermite : public FEPhysicalElement<t_name, t_npts>
{
public:
    FEHermite ();
    ~FEHermite () override;

private:
    FEHermite (const FEHermite &) = delete;
    FEHermite & operator= (const FEHermite &) = delete;
};

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
FEHermite<t_name, t_npts, t_nddl>::FEHermite () : FEPhysicalElement<t_name, t_npts> ()
{
    this->order            = 0;
    this->dim              = 0;
    this->namePhysicalBase = "Emp";
}

template <PHYSICAL_CELL_TYPE t_name, int t_npts, int t_nddl>
FEHermite<t_name, t_npts, t_nddl>::~FEHermite ()
{
}

#endif /* SRC_SOLVER_FESTRUCT_LISTOFELEMENTS_FEHERMITE_FEHERMITE_HPP */
