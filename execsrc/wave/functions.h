#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <Soliton>

class FunAnalytic : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunAnalytic)

    real_t ToReal(Point *p, real_t t, Cell *cell) override;
};

class FunAnalytic_1 : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunAnalytic_1)

    real_t ToReal(Point *p, real_t t, Cell *cell) override;
};

class FunAnalytic_2 : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunAnalytic_2)

    real_t ToReal(Point *p, real_t t, Cell *cell) override;
};

class FunNeumann : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunNeumann)

    Point To3DPoint(Point *p, real_t t, Cell *cell) override;
};

class FunDirichlet : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunDirichlet)

    real_t ToReal(Point *p, real_t t, Cell *cell) override;
};

class FunSecondMember : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunSecondMember)

    real_t ToReal(Point *p, real_t t, Cell *cell) override;
};

class FunJumpZeta : public SolitonFunctor
{
public:
    SOLITON_FUNCTOR_DEF(FunJumpZeta)

    real_t ToReal(Point *p, real_t t, Cell *cell) override;
};

extern FunAnalytic      fun_analytic;
extern FunAnalytic_1    fun_analytic_1;
extern FunAnalytic_2    fun_analytic_2;
extern FunNeumann       fun_neumann;
extern FunDirichlet     fun_dirichlet;
extern FunSecondMember  fun_secondmember;
extern FunJumpZeta      fun_jumpzeta;
extern real_t           coeffA;
extern real_t           coeffT;
extern real_t           coeffLambda;

#endif // FUNCTIONS_H
