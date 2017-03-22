#ifndef __BASIS_H
#define __BASIS_H

#include <boost/dynamic_bitset.hpp>
#include <cmath>

#include "Environment.h"

class Basis
{
  private:
    unsigned int l_, n_;
    LLInt factorial_(LLInt n);
    LLInt first_int_();

  public:
    LLInt basis_size;
    PetscInt basis_local;
    PetscInt basis_start;
    PetscInt nlocal;
    PetscInt start;
    PetscInt end;
    LLInt *int_basis;
    Basis(const Environment &env);
    ~Basis();
    Basis(const Basis &rhs);
    Basis &operator=(const Basis &rhs);
    void construct_int_basis();
    void print_basis(const Environment &env);
    void construct_bit_basis(boost::dynamic_bitset<> *bit_basis);
};
#endif
