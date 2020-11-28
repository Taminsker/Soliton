#ifndef ROUTINES_H
#define ROUTINES_H

#include <Core/core.h>
#include <Algorithms/algorithms.h>
#include <Solver/solver.h>
#include <IO/io.h>

namespace SolitonRoutines {
void PrintMacros ();
void PrintExecutableName (std::string executablename);

void Generation (MeshStorage* store, InputDatStruct* inputdatfile);

void ComputeGradGrad (Mesh* mesh, std::vector<Triplet_eig>* listTriplets);

void ComputeSecondMember (Mesh* mesh, PlainVector_eig *vec, real_t (*f)(Point, real_t));


void ComputeDirichletNoSBM (Mesh* mesh, std::vector<Triplet_eig>* listTriplets, PlainVector_eig *F, real_t (*f)(Point, real_t), PHYS tagToApply, real_t hPerp = -1, real_t t = 0);

void ComputeDirichletSBM (Mesh* mesh, Mesh* object, std::vector<Triplet_eig>* listTriplets, PlainVector_eig *F, real_t (*f)(Point, real_t), real_t hPerp = -1, real_t t = 0);

void ComputeNeumannNoSBM (Mesh* mesh, std::vector<Triplet_eig>* listTriplets, PlainVector_eig *F, real_t (*f)(Point, real_t), PHYS tagToApply, real_t hPerp = 1., real_t t = 0);


void PostProcess (MeshStorage *store, std::vector<SolitonFEContainer*> *containers, std::vector<PlainVector_eig*> *sols, std::function<real_t(Point, real_t)> fun);

}
#endif // ROUTINES_H
