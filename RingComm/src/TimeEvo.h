#ifndef __TIME_EVO_H
#define __TIME_EVO_H

#include "Environment.h"
#include "Basis.h"
#include "SparseOp.h"

class TimeEvo
{
  private:
    unsigned int l_;
    unsigned int n_;
    PetscMPIInt mpirank_;
    PetscMPIInt mpisize_;

  public:
    TimeEvo(const Environment &env);
    ~TimeEvo();
    void loschmidt_echo(const SparseOp &sparse, 
                        const unsigned int &iterations, 
                        const double *times, 
                        const double &tol, 
                        const int &maxits, 
                        const Vec &v, 
                        const Mat &ham_mat);
};
#endif
