#ifndef SRC_ROUTINES_ROUTINES_HPP
#define SRC_ROUTINES_ROUTINES_HPP

#include "../Algorithms/algorithms.hpp"
#include "../Core/core.hpp"
#include "../IO/io.hpp"
#include "../Solver/solver.hpp"

namespace SolitonRoutines
{
void PrintMacros ();
void PrintExecutableName (std::string executablename);

void Generation (MeshStorage * store, InputDatStruct * inputdatfile);

void ComputeGradGrad (Mesh * mesh, std::vector<Triplet> * listTriplets);

void ComputeSecondMember (Mesh * mesh, DenseVector * vec, real_t (*f) (Point, real_t));

void ComputeDirichletNoSBM (Mesh * mesh, std::vector<Triplet> * listTriplets, DenseVector * F, real_t (*f) (Point, real_t), PHYS tagToApply, real_t hPerp = -1, real_t t = 0);

void ComputeDirichletSBM (Mesh * mesh, Mesh * object, std::vector<Triplet> * listTriplets, DenseVector * F, real_t (*f) (Point, real_t), real_t hPerp = -1, real_t t = 0);

void ComputeNeumannNoSBM (Mesh * mesh, std::vector<Triplet> * listTriplets, DenseVector * F, real_t (*f) (Point, real_t), PHYS tagToApply, real_t hPerp = 1., real_t t = 0);

void PostProcess (MeshStorage * store, std::vector<SolitonFEContainer *> * containers, std::vector<DenseVector *> * sols, std::function<real_t (Point, real_t)> fun);

}  // namespace SolitonRoutines
#endif /* SRC_ROUTINES_ROUTINES_HPP */
