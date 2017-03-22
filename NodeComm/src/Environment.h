#ifndef __ENVIRONMENT_H
#define __ENVIRONMENT_H

#include <iostream>
#include <stdint.h>

#include <petscsys.h>
#include <slepcmfn.h>
#include <slepceps.h>

typedef unsigned long long ULLInt;
typedef PetscInt LLInt;

class Environment
{
  private:

  public:
    unsigned int l, n;
    PetscMPIInt mpirank;
    PetscMPIInt mpisize;
    PetscMPIInt node_rank;
    PetscMPIInt node_size;
    MPI_Comm node_comm;
    Environment(int argc, char **argv, unsigned int l, unsigned int n);
    ~Environment();
    LLInt basis_size() const;
    void distribution(PetscInt b_size, 
                      PetscInt &nlocal, 
                      PetscInt &start, 
                      PetscInt &end) const;
};
#endif
