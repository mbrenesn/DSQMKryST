#ifndef __SPARSEOP_H
#define __SPARSEOP_H

#include <boost/math/constants/constants.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <stdexcept>
#include <ctime>

#include "Basis.h"
#include "Environment.h"

class SparseOp
{
  private:
    unsigned int l_, n_;
    PetscMPIInt mpirank_;
    PetscMPIInt mpisize_;
    LLInt basis_size_;
    PetscInt nlocal_;
    PetscInt start_;
    PetscInt end_;
    LLInt mod_(LLInt a, LLInt b);
    void gather_nonlocal_values_(LLInt *start_inds);
    void determine_allocation_details_(LLInt *int_basis, 
                                       std::vector<LLInt> &cont,
                                       std::vector<LLInt> &st, 
                                       PetscInt *diag, 
                                       PetscInt *off);

  public:
    SparseOp(const Environment &env, const Basis &basis);
    ~SparseOp();
    inline LLInt binary_to_int(boost::dynamic_bitset<> bs);
    inline LLInt binsearch(const LLInt *array, LLInt len, LLInt value);
    LLInt get_neel_index(const Basis &bas);
    LLInt get_random_index(const Basis &bas, bool wtime, bool verbose);
    void construct_AA_hamiltonian(Mat &ham_mat, 
                                  LLInt *int_basis, 
                                  double V,
                                  double t, 
                                  double h,
                                  double beta);
};
#endif
