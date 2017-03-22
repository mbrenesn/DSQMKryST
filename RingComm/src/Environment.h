#ifndef __ENVIRONMENT_H
#define __ENVIRONMENT_H

#include <iostream>

#include <petscsys.h>
#include <slepcmfn.h>
#include <slepceps.h>

typedef unsigned long long int ULLInt;
typedef PetscInt LLInt;

class Environment
{
  private:

  public:
    unsigned int l, n;
    PetscMPIInt mpirank;
    PetscMPIInt mpisize;
    Environment(int argc, char **argv, unsigned int l, unsigned int n);
    ~Environment();
    LLInt basis_size() const;
    void distribution(PetscInt b_size, 
                      PetscInt &nlocal, 
                      PetscInt &start, 
                      PetscInt &end) const;
};
#endif
