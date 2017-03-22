#include "SparseOp.h"

SparseOp::SparseOp(const Environment &env, const Basis &basis)
{
  l_ = env.l;
  n_ = env.n;
  mpirank_ = env.mpirank;
  mpisize_ = env.mpisize;
  nlocal_ = basis.nlocal;
  start_ = basis.start;
  end_ = basis.end;
  basis_size_ = basis.basis_size;
}

SparseOp::~SparseOp()
{}

LLInt SparseOp::mod_(LLInt a, LLInt b)
{
  return (a % b + b) % b;
}

inline
LLInt SparseOp::binary_to_int(boost::dynamic_bitset<> bs)
{
  LLInt integer = 0;

  for(unsigned int i = 0; i < l_; ++i){
    if(bs[i] == 1){
      integer += 1ULL << i;
    }
  }

  return integer;
}

/*******************************************************************************/
// Binary search: Divide and conquer. For the construction of the Hamiltonian
// matrix instead of looking through all the elements of the int basis a
// binary search will perform better for large systems
/*******************************************************************************/
inline
LLInt SparseOp::binsearch(const LLInt *array, LLInt len, LLInt value)
{
  if(len == 0) return -1;
  LLInt mid = len / 2;

  if(array[mid] == value) 
    return mid;
  else if(array[mid] < value){
    LLInt result = binsearch(array + mid + 1, len - (mid + 1), value);
    if(result == -1) 
      return -1;
    else
      return result + mid + 1;
  }
  else
    return binsearch(array, mid, value);
}

/*******************************************************************************/
// Initial Neel state index 
/*******************************************************************************/
LLInt SparseOp::get_neel_index(const Basis &bas)
{
  LLInt index;
  if(l_ / 2 != n_){
    std::cerr << "Not implemented!" << std::endl;
    std::cerr << "Neel state has only been implemented for half-filled systems" << std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  boost::dynamic_bitset<> neel(l_, 1);
  for(unsigned int site = 0; site < l_; site += 2){
    neel.set(site);
  }

  LLInt neel_int = binary_to_int(neel);
  index = binsearch(bas.int_basis, bas.basis_local, neel_int);
  if(index != -1) index += start_;

  return index;
}

/*******************************************************************************/
// Pick initial index out of the basis randomly 
/*******************************************************************************/
LLInt SparseOp::get_random_index(const Basis &bas, bool wtime, bool verbose)
{
  LLInt pick_ind;
  boost::random::mt19937 gen;

  if(wtime) gen.seed(static_cast<LLInt>(std::time(0)));

  if(mpirank_ == 0){
    boost::random::uniform_int_distribution<LLInt> dist(0, bas.basis_size - 1);
    pick_ind = dist(gen);
  }
  MPI_Bcast(&pick_ind, 1, MPI_LONG_LONG_INT, 0, PETSC_COMM_WORLD);

  bool check = false;
  if(pick_ind >= start_ && pick_ind < end_) check = true;
    
  if(verbose && check){
    std::cout << "Initial state randomly chosen: " << bas.int_basis[pick_ind] << std::endl;
    std::cout << "With binary representation: " << std::endl;
    boost::dynamic_bitset<> bs(l_, bas.int_basis[pick_ind - start_]);
    std::cout << bs << std::endl;
  }

  return pick_ind;
}

/*******************************************************************************/
// Collects the nonlocal start index (global indices) of every processor, to be used
// during the construction of the matrix
/*******************************************************************************/
void SparseOp::gather_nonlocal_values_(LLInt *start_inds)
{
  MPI_Allgather(&start_, 1, MPI_LONG_LONG, start_inds, 1, MPI_LONG_LONG, 
    PETSC_COMM_WORLD);
}

/*******************************************************************************/
// Determines the sparsity pattern to allocate memory only for the non-zero 
// entries of the matrix
/*******************************************************************************/
void SparseOp::determine_allocation_details_(LLInt *int_basis, 
                                             std::vector<LLInt> &cont, 
                                             std::vector<LLInt> &st, 
                                             PetscInt *diag, 
                                             PetscInt *off)
{
  for(PetscInt i = 0; i < nlocal_; ++i) diag[i] = 1;

  for(PetscInt state = start_; state < end_; ++state){

    boost::dynamic_bitset<> bs(l_, int_basis[state - start_]);

    // Loop over all sites of the bit representation
    for(unsigned int site = 0; site < l_; ++site){
      // A copy to avoid modifying the original basis
      boost::dynamic_bitset<> bitset = bs;
      
      // Case 1: There's a particle in this site
      if(bitset[site] == 1){
        int next_site1 = (site + 1) % l_;

        // If there's a particle in next site, do nothing
        if(bitset[next_site1] == 1){
          continue;
        }
        // Otherwise do a swap
        else{
          bitset[next_site1] = 1;
          bitset[site]       = 0;

          LLInt new_int1 = binary_to_int(bitset);
          // Loop over all states and look for a match
          LLInt match_ind1 = binsearch(int_basis, nlocal_, new_int1);
          if(match_ind1 == -1){
            cont.push_back(new_int1);
            st.push_back(state);
            continue;
          }
          else{
            match_ind1 += start_;
          }

          if(match_ind1 < end_ && match_ind1 >= start_) diag[state - start_]++;
          else off[state - start_]++; 
        }
      }
      // Case 2: There's no particle in this site
      else{
        int next_site0 = (site + 1) % l_;

        // If there's a particle in the next site, a swap can occur
        if(bitset[next_site0] == 1){
          bitset[next_site0] = 0;
          bitset[site]       = 1;

          LLInt new_int0 = binary_to_int(bitset);
          // Loop over all states and look for a match
          LLInt match_ind0 = binsearch(int_basis, nlocal_, new_int0);
          if(match_ind0 == -1){
            cont.push_back(new_int0);
            st.push_back(state);
            continue;
          }
          else{
            match_ind0 += start_;
          }
          
          if(match_ind0 < end_ && match_ind0 >= start_) diag[state - start_]++;
          else off[state - start_]++;
        }
        // Otherwise do nothing
        else{
          continue;
        }
      }    
    }
  }

  // Collective communication of global indices
  LLInt *start_inds = new LLInt[mpisize_];

  gather_nonlocal_values_(start_inds);

  // Proc 0 is always gonna have the larger section of the distribution (when rest is present)
  // so let's use this value as the size of the basis_help buffer
  LLInt basis_help_size;
  if(mpirank_ == 0) basis_help_size = nlocal_;

  MPI_Bcast(&basis_help_size, 1, MPI_LONG_LONG, 0, PETSC_COMM_WORLD);

  // Create basis_help buffers and initialize them to zero
  LLInt *basis_help = new LLInt[basis_help_size];
  for(LLInt i = 0; i < basis_help_size; ++i) basis_help[i] = 0;

  // At the beginning basis_help is just int_basis, with the remaining values set to zero
  // It's important that the array remains sorted for the binary lookup
  for(LLInt i = 0; i < nlocal_; ++i)
    basis_help[i + (basis_help_size - nlocal_)] = int_basis[i];

  // Main communication procedure. A ring exchange of the int_basis using basis_help memory
  // buffer, after a ring exchange occurs each processor looks for the missing indices of
  // the Hamiltonian and replaces the values of cont buffer
  PetscMPIInt next = (mpirank_ + 1) % mpisize_;
  PetscMPIInt prec = (mpirank_ + mpisize_ - 1) % mpisize_;

  LLInt cont_size = cont.size();

  for(PetscMPIInt exc = 0; exc < mpisize_ - 1; ++exc){
 
    MPI_Sendrecv_replace(&basis_help[0], basis_help_size, MPI_LONG_LONG, next, 0,
      prec, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

    PetscMPIInt source = mod_((prec - exc), mpisize_);
    for(LLInt i = 0; i < cont_size; ++i){
      if(cont[i] > 0){
        LLInt m_ind = binsearch(basis_help, basis_help_size, cont[i]);

        if(m_ind != -1){
          if(basis_help[0] == 0) m_ind = m_ind - 1;
          cont[i] = -1ULL * (m_ind + start_inds[source]);
        }
      }
    }
  }
  
  // Flip the signs
  for(PetscInt i = 0; i < cont_size; ++i) cont[i] = -1ULL * cont[i];
  
  //std::cout << "Cont from Proc " << mpirank_ << std::endl;
  //for(LLInt i = 0; i < cont_size; ++i) std::cout << cont[i] << std::endl;

  // Now cont contains the missing indices
  for(LLInt in = 0; in < cont_size; ++in){
    LLInt st_c = st[in];
    if(cont[in] < end_ && cont[in] >= start_) diag[st_c - start_]++;
    else off[st_c - start_]++;
  }

  delete [] basis_help;
  delete [] start_inds;
}

/*******************************************************************************/
// Computes the Hamiltonian matrix given by means of the integer basis
/*******************************************************************************/
void SparseOp::construct_AA_hamiltonian(Mat &ham_mat, 
                                        LLInt *int_basis, 
                                        double V, 
                                        double t, 
                                        double h,
                                        double beta)
{
  // Preallocation. For this we need a hint on how many non-zero entries the matrix will
  // have in the diagonal submatrix and the offdiagonal submatrices for each process

  // Allocating memory only for the non-zero entries of the matrix
  PetscInt *d_nnz, *o_nnz;
  PetscCalloc1(nlocal_, &d_nnz);
  PetscCalloc1(nlocal_, &o_nnz);

  std::vector<LLInt> cont;
  cont.reserve(basis_size_ / l_);
  std::vector<LLInt> st;
  st.reserve(basis_size_ / l_);
 
  determine_allocation_details_(int_basis, cont, st, d_nnz, o_nnz);

  // Preallocation step
  MatMPIAIJSetPreallocation(ham_mat, 0, d_nnz, 0, o_nnz);

  PetscFree(d_nnz);
  PetscFree(o_nnz);

  // Hamiltonian matrix construction
  PetscScalar ti = t;
  PetscScalar Vi = V;
  const double pi = boost::math::constants::pi<double>();

  // Grab 1 of the states and turn it into bit representation
  for(PetscInt state = start_; state < end_; ++state){ 
    
    boost::dynamic_bitset<> bs(l_, int_basis[state - start_]);

    // Loop over all sites of the bit representation
    for(unsigned int site = 0; site < l_; ++site){
      // A copy to avoid modifying the original basis
      boost::dynamic_bitset<> bitset = bs;
  
      // Case 1: There's a particle in this site
      if(bitset[site] == 1){
        
        PetscScalar osc_term = h * cos(2 * pi * beta * site);
        MatSetValues(ham_mat, 1, &state, 1, &state, &osc_term, ADD_VALUES);
        
        int next_site1 = (site + 1) % l_;

        // If there's a particle in next site, do nothing
        if(bitset[next_site1] == 1){
          // Accumulate 'V' terms
          MatSetValues(ham_mat, 1, &state, 1, &state, &Vi, ADD_VALUES);
          continue;
        }
        // Otherwise do a swap
        else{
          bitset[next_site1] = 1;
          bitset[site]       = 0;

          LLInt new_int1 = binary_to_int(bitset);
          // Loop over all states and look for a match
          LLInt match_ind1 = binsearch(int_basis, nlocal_, new_int1);
          if(match_ind1 == -1){
            continue;
          }
          else{
            match_ind1 += start_;
          }
          
          if(match_ind1 == -1){
            std::cerr << "Error in the binary search within the Ham mat construction" 
              << std::endl;
            MPI_Abort(PETSC_COMM_WORLD, 1);
          } 
            
          MatSetValues(ham_mat, 1, &match_ind1, 1, &state, &ti, ADD_VALUES);
        }
      }
      // Case 2: There's no particle in this site
      else{
        int next_site0 = (site + 1) % l_;

        // If there's a particle in the next site, a swap can occur
        if(bitset[next_site0] == 1){
          bitset[next_site0] = 0;
          bitset[site]       = 1;

          LLInt new_int0 = binary_to_int(bitset);
          // Loop over all states and look for a match
          LLInt match_ind0 = binsearch(int_basis, nlocal_, new_int0);
          if(match_ind0 == -1){
            continue;
          }
          else{
            match_ind0 += start_;
          }
          
          if(match_ind0 == -1){
            std::cerr << "Error in the binary search within the Ham mat construction" 
              << std::endl;
            MPI_Abort(PETSC_COMM_WORLD, 1);
          } 
              
          MatSetValues(ham_mat, 1, &match_ind0, 1, &state, &ti, ADD_VALUES);
        }
        // Otherwise do nothing
        else{
          continue;
        }
      }    
    }
  }

  // Cont already contains the missing indices
  for(ULLInt in = 0; in < cont.size(); ++in){
    LLInt st_c = st[in];
    LLInt cont_c = cont[in];
    MatSetValues(ham_mat, 1, &cont_c, 1, &st_c, &ti, ADD_VALUES);
  }

  MatAssemblyBegin(ham_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(ham_mat, MAT_FINAL_ASSEMBLY);

  MatSetOption(ham_mat, MAT_SYMMETRIC, PETSC_TRUE);
}
