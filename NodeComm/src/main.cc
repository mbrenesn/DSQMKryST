#include <iostream>

#include "Environment.h"
#include "Basis.h"
#include "SparseOp.h"
#include "TimeEvo.h"

#include <petsctime.h>

int main(int argc, char **argv)
{
  unsigned int l = 55;
  unsigned int n = 1;
  double V = 1.0;
  double t = 0.5;
  double h = 1.0;
  double beta = 34.0 / 55.0;

  // Establish the environment
  Environment env(argc, argv, l, n);

  // Establish the basis environment, by pointer, to call an early destructor and reclaim
  // basis memory
  Basis *basis = new Basis(env);

  // Construct basis
  basis->construct_int_basis();
  //basis->print_basis(env);

  // Establish Hamiltonian operator environment
  SparseOp aubry(env, *basis);

  // Create a matrix object that will be used as Hamiltonian
  Mat aa_mat;
  MatCreate(PETSC_COMM_WORLD, &aa_mat);
  MatSetSizes(aa_mat, basis->nlocal, basis->nlocal, basis->basis_size, basis->basis_size);
  MatSetType(aa_mat, MATMPIAIJ);
  
  // Construct the Hamiltonian matrix
  aubry.construct_AA_hamiltonian(aa_mat, 
                                 basis->int_basis, 
                                 V, 
                                 t, 
                                 h,
                                 beta);

  // Neel index before destroying basis
  //LLInt neel_index = aubry.get_neel_index(*basis);

  // Random index before destroying basis
  LLInt random_index = aubry.get_random_index(*basis, true, true);

  // Call basis destructor to reclaim memory
  delete basis;

  // Initial state with the same parallel layout as the matrix
  Vec v;
  MatCreateVecs(aa_mat, NULL, &v);
  VecZeroEntries(v);
  //VecSetValue(v, neel_index, 1.0, INSERT_VALUES);
  VecSetValue(v, random_index, 1.0, INSERT_VALUES);
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);

  /*** Time evolultion ***/
  PetscLogDouble kryt1, kryt2;
  PetscTime(&kryt1);

  TimeEvo te(env);

  double tol = 1.0e-7;
  int maxits = 100000;

  const int iterations = 36;
  double times[iterations + 1] 
      = {0.0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,
         0.6,0.7,0.8,0.9,1.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,
         45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0};

  te.loschmidt_echo(aubry, iterations, times, tol, maxits, v, aa_mat);
 
  PetscTime(&kryt2);
  /*** End of time evolution ***/

  VecDestroy(&v);
  MatDestroy(&aa_mat);

  return 0;
}
