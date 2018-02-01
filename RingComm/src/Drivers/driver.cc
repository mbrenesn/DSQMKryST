/** @addtogroup RingComm */
/** @file */
#include <iostream>

#include "../Environment/Environment.h"
#include "../Utils/Utils.h"
#include "../Basis/Basis.h"
#include "../Operators/SparseOp.h"
#include "../InitialState/InitialState.h"
#include "../TimeEvo/KrylovEvo.h"

int main(int argc, char **argv)
{
  unsigned int l = 55;
  unsigned int n = 1;
  double V = 1.0;
  double t = 0.5;
  double h = 1.0;
  double beta = 34.0 / 55.0;
  double initial_time = 0.0;
  double final_time = 10.0;

  // Establish the environment
  EnvironmentRC env(argc, argv, l, n);

  PetscMPIInt mpirank = env.mpirank;
  PetscMPIInt mpisize = env.mpisize;

  // Establish the basis environment, by pointer, to call an early destructor and reclaim
  // basis memory
  BasisRC *basis = new BasisRC(env);

  // Construct basis
  basis->construct_int_basis();
  //basis->print_basis(env);

  // Establish the Hamiltonian operator environment
  SparseOpRC aubry(env, *basis); 

  // Construct the Hamiltonian matrix
  aubry.construct_AA_hamiltonian(basis->int_basis,
                                 V,
                                 t, 
                                 h,
                                 beta);

  // Create an initial state before deleting the basis
  InitialStateRC init(env, *basis);
  init.random_initial_state(basis->int_basis, false, true);

  delete basis;

  // Time Evo, the Krylov subspace method based on Arnoldi decomposition is invoked here
  double tol = 1.0e-7;
  int maxits = 1000000;

  KrylovEvoRC te(aubry.HamMat, tol, maxits);

  // Initial value
  PetscScalar l_echo;
  PetscReal ld;
  Vec t0_vec;
  VecDuplicate(init.InitialVec, &t0_vec);
  VecCopy(init.InitialVec, t0_vec);
  VecDot(t0_vec, init.InitialVec, &l_echo);

  ld = (PetscRealPart(l_echo) * PetscRealPart(l_echo)) + 
      (PetscImaginaryPart(l_echo) * PetscImaginaryPart(l_echo));
  if(mpirank == 0){
    std::cout << "Time"  << "\t" << "Loschmidt echo" << std::endl;
    std::cout << "0.0" << "\t" << ld << std::endl;
  }

  // Time evo
  te.krylov_evo(final_time, initial_time, init.InitialVec);
  VecDot(t0_vec, init.InitialVec, &l_echo);
  ld = (PetscRealPart(l_echo) * PetscRealPart(l_echo)) + 
      (PetscImaginaryPart(l_echo) * PetscImaginaryPart(l_echo));
  if(mpirank == 0){
    std::cout << final_time << "\t" << ld << std::endl;
  }

  VecDestroy(&t0_vec);
  return 0;
}
