#ifndef ROUTINES_H
#define ROUTINES_H

#include <Core/core.h>
#include <Algorithms/algorithms.h>
#include <Solver/solver.h>
#include <IO/io.h>

namespace SolitonRoutines {
void PrintMacros ();
void PrintExecutableName (std::string executablename);

void Generation (Sto4Sol* store, InputDatStruct* inputdatfile);

void ComputeGradGrad (Mesh* mesh, std::vector<Triplet>* listTriplets);

void ComputeSecondMember (Mesh* mesh, PlainVector* vec, double (*f)(Point, double));


void ComputeDirichletNoSBM (Mesh* mesh, std::vector<Triplet>* listTriplets, PlainVector* F, double (*f)(Point, double), PHYS tagToApply, double hPerp = -1, double t = 0);

void ComputeDirichletSBM (Mesh* mesh, Mesh* object, std::vector<Triplet>* listTriplets, PlainVector* F, double (*f)(Point, double), double hPerp = -1, double t = 0);

void ComputeNeumannNoSBM (Mesh* mesh, std::vector<Triplet>* listTriplets, PlainVector* F, double (*f)(Point, double), PHYS tagToApply, double hPerp = 1., double t = 0);


void PostProcess (Sto4Sol *store, SuperSolver* solver, std::vector<PlainVector*>* sols, std::function<double(Point, double)> fun);

}
#endif // ROUTINES_H
