#include "Environment.h"

EnvironmentRC::EnvironmentRC(int argc, 
                           char **argv, 
                           unsigned int l, 
                           unsigned int n)
: l(l), n(n)
{
  SlepcInitialize(&argc, &argv, NULL, NULL);

  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank);
}

EnvironmentRC::~EnvironmentRC()
{
  SlepcFinalize();
}

/*******************************************************************************/
// Computes the size of the Hilbert space by means of using all possible com
// binations, this would mean using the factorial function as:
// return factorial_(l_) / ( factorial_(n_) * factorial_(l_ - n_) );
// The problem with this is that datatype overflowing can occur even for uns
// igned long long integers. Instead we can use the following expression to
// compute the size of the system.
/*******************************************************************************/
LLInt EnvironmentRC::basis_size() const 
{
  double size = 1.0;
  for(LLInt i = 1; i <= (l - n); ++i){
    size *= (static_cast<double> (i + n) / static_cast<double> (i));  
  }

  return floor(size + 0.5);
}

/*******************************************************************************/
// Determines the distribution among processes without relying on PETSc routines
// This is done to avoid allocation of PETSc objects before they are required,
// therefore saving memory
/*******************************************************************************/
void EnvironmentRC::distribution(PetscInt b_size, 
                                 PetscInt &nlocal, 
                                 PetscInt &start, 
                                 PetscInt &end) const
{
  nlocal = b_size / mpisize;
  PetscInt rest = b_size % mpisize;

  if(rest && (mpirank < rest)) nlocal++;

  start = mpirank * nlocal;
  if(rest && (mpirank >= rest)) start += rest;

  end = start + nlocal;
}
