#ifndef SOLITONFEITEM_H
#define SOLITONFEITEM_H

#include "../SolitonFEItemBase/solitonfeitembase.h"

template<ITEM_T t_type>
class SolitonFEItem : public SolitonFEItemBase
{
public:
    SolitonFEItem (SolitonFEContainer *parent, PHYS tag = PHYS::DEFAULT) :
        SolitonFEItemBase (parent, tag)
    {}

    SolitonFEItem (const SolitonFEItem& tocopy) :
        SolitonFEItemBase (tocopy)
    {}

    ~SolitonFEItem () override {}

    void Compute (real_t time) override;

    SOLITON_INLINE
    void SetAdditionalName (std::string name) override
    {
        this->m_name = to_string (t_type) + name;
        return;
    }

    SOLITON_INLINE
    ITEM_T GetType() const override
    {
        return t_type;
    }
};


#define EXPLICIT_INST(X) \
    template<>\
    void\
    SolitonFEItem<ITEM_T::EMPTY>::Compute (real_t time)

EXPLICIT_INST(EMPTY);
EXPLICIT_INST(PHI_PHI);
EXPLICIT_INST(DAMPING_PHI_PHI);
EXPLICIT_INST(GRAD_GRAD);
EXPLICIT_INST(SECOND_MEMBER);
EXPLICIT_INST(DIRICHLET_BOUNDARY);
EXPLICIT_INST(NEUMANN_BOUNDARY);
EXPLICIT_INST(DIRICHLET_BOUNDARY_SBM);
EXPLICIT_INST(NEUMANN_BOUNDARY_SBM);


#undef EXPLICIT_INST

template<>
void
SolitonFEItem<ITEM_T::EMPTY>::Compute (real_t time);

#endif // SUPERFEITEM_H
