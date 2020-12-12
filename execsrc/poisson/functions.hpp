#ifndef EXECSRC_POISSON_FUNCTIONS_HPP
#define EXECSRC_POISSON_FUNCTIONS_HPP

#include <Soliton>

class FunAnalytic : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF (FunAnalytic)

    real_t ToReal (Point * p, real_t t, Cell * cell) override;
};

class FunNeumann : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF (FunNeumann)

    Point To3DPoint (Point * p, real_t t, Cell * cell) override;
};

class FunDirichlet : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF (FunDirichlet)

    real_t ToReal (Point * p, real_t t, Cell * cell) override;
};

class FunSecondMember : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF (FunSecondMember)

    real_t ToReal (Point * p, real_t t, Cell * cell) override;
};

extern FunAnalytic     fun_analytic;
extern FunNeumann      fun_neumann;
extern FunDirichlet    fun_dirichlet;
extern FunSecondMember fun_secondmember;

#endif /* EXECSRC_POISSON_FUNCTIONS_HPP */
