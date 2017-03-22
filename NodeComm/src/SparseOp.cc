#include "SparseOp.h"

SparseOp::SparseOp(const Environment &env, const Basis &basis)
{
  l_ = env.l;
  n_ = env.n;
  mpirank_ = env.mpirank;
  mpisize_ = env.mpisize;
  node_rank_ = env.node_rank;
  node_size_ = env.node_size;
  nlocal_ = basis.nlocal;
  start_ = basis.start;
  end_ = basis.end;
  basis_size_ = basis.basis_size;
  MPI_Comm_dup(env.node_comm, &node_comm_);
}

SparseOp::~SparseOp()
{
  MPI_Comm_free(&node_comm_);
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
    std::cerr << "Neel state has only been implemented for half-filled systems" 
      << std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  boost::dynamic_bitset<> neel(l_, 1);
  for(unsigned int site = 0; site < l_; site += 2){
    neel.set(site);
  }

  if(mpirank_ == 0){
    LLInt neel_int = binary_to_int(neel);
    index = binsearch(bas.int_basis, bas.basis_local, neel_int);
  }

  MPI_Bcast(&index, 1, MPI_LONG_LONG_INT, 0, PETSC_COMM_WORLD);

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

    if(verbose){
      std::cout << "Initial state randomly chosen: " << bas.int_basis[pick_ind] 
        << std::endl;
      std::cout << "With binary representation: " << std::endl;
      boost::dynamic_bitset<> bs(l_, bas.int_basis[pick_ind]);
      std::cout << bs << std::endl;
    }
  }
  
  MPI_Bcast(&pick_ind, 1, MPI_LONG_LONG_INT, 0, PETSC_COMM_WORLD);

  return pick_ind;
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

    PetscInt basis_ind;
    node_rank_ ? basis_ind = state - start_ : basis_ind = state;

    boost::dynamic_bitset<> bs(l_, int_basis[basis_ind]);

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
          LLInt match_ind1;
          if(node_rank_){
            match_ind1 = binsearch(int_basis, nlocal_, new_int1); 
            if(match_ind1 == -1){
              cont.push_back(new_int1);
              st.push_back(state);
              continue;
            }
            else{
              match_ind1 += start_;
            }
          }
          else{
            match_ind1 = binsearch(int_basis, basis_size_, new_int1);
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
          LLInt match_ind0;
          if(node_rank_){
            match_ind0 = binsearch(int_basis, nlocal_, new_int0); 
            if(match_ind0 == -1){
              cont.push_back(new_int0);
              st.push_back(state);
              continue;
            }
            else{
              match_ind0 += start_;
            }
          }
          else{
            match_ind0 = binsearch(int_basis, basis_size_, new_int0);
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

  LLInt *recv_sizes = NULL;
  if(node_rank_ == 0) recv_sizes = new LLInt[node_size_ - 1];
  // Communication to rank 0 of every node to find size of buffers
  if(node_rank_){
    LLInt cont_size = cont.size();
    MPI_Send(&cont_size, 1, MPI_LONG_LONG, 0, node_rank_, node_comm_);
  }
  else{
    for(PetscMPIInt i = 1; i < node_size_; ++i)
      MPI_Recv(&recv_sizes[i - 1], 1, MPI_LONG_LONG, i, MPI_ANY_TAG, node_comm_, 
       MPI_STATUS_IGNORE);
  }

  // Communication to rank 0 of each node to find missing indices
  if(node_rank_){
    MPI_Send(&cont[0], cont.size(), MPI_LONG_LONG, 0, node_rank_, node_comm_);
    MPI_Recv(&cont[0], cont.size(), MPI_LONG_LONG, 0, 0, node_comm_, MPI_STATUS_IGNORE);
  }
  else{
    for(PetscMPIInt i = 1; i < node_size_; ++i){
      MPI_Status stat;
      LLInt rsize = recv_sizes[i - 1];
      cont.resize(rsize);
      MPI_Recv(&cont[0], recv_sizes[i - 1], MPI_LONG_LONG, i, i, node_comm_,
        &stat);
    
      for(LLInt ii = 0; ii < rsize; ++ii){ 
        LLInt m_ind = binsearch(int_basis, basis_size_, cont[ii]);
        cont[ii] = m_ind;
      }

      MPI_Send(&cont[0], recv_sizes[i - 1], MPI_LONG_LONG, stat.MPI_SOURCE, 0, node_comm_);
      cont.erase(cont.begin(), cont.end());
    }
  }

  // Now cont contains the missing indices
  if(node_rank_){
    for(ULLInt in = 0; in < cont.size(); ++in){
      LLInt st_c = st[in];
      if(cont[in] < end_ && cont[in] >= start_) diag[st_c - start_]++;
      else off[st_c - start_]++;
    }
  }
  
  delete [] recv_sizes;
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
    
    PetscInt basis_ind;
    node_rank_ ? basis_ind = state - start_ : basis_ind = state;
    
    boost::dynamic_bitset<> bs(l_, int_basis[basis_ind]);

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
          LLInt match_ind1;
          if(node_rank_){
            match_ind1 = binsearch(int_basis, nlocal_, new_int1); 
            if(match_ind1 == -1){
              continue;
            }
            else{
              match_ind1 += start_;
              MatSetValues(ham_mat, 1, &match_ind1, 1, &state, &ti, ADD_VALUES);
            }
          }
          else{
            match_ind1 = binsearch(int_basis, basis_size_, new_int1);
            MatSetValues(ham_mat, 1, &match_ind1, 1, &state, &ti, ADD_VALUES);
          }
          
          if(match_ind1 == -1){
            std::cerr << "Error in the binary search within the Ham mat construction" 
              << std::endl;
            MPI_Abort(PETSC_COMM_WORLD, 1);
          } 
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
          LLInt match_ind0;
          if(node_rank_){
            match_ind0 = binsearch(int_basis, nlocal_, new_int0); 
            if(match_ind0 == -1){
              continue;
            }
            else{
              match_ind0 += start_;
              MatSetValues(ham_mat, 1, &match_ind0, 1, &state, &ti, ADD_VALUES);
            }
          }
          else{
            match_ind0 = binsearch(int_basis, basis_size_, new_int0);
            MatSetValues(ham_mat, 1, &match_ind0, 1, &state, &ti, ADD_VALUES);
          }
          
          if(match_ind0 == -1){
            std::cerr << "Error in the binary search within the Ham mat construction" 
              << std::endl;
            MPI_Abort(PETSC_COMM_WORLD, 1);
          } 
        }
        // Otherwise do nothing
        else{
          continue;
        }
      }    
    }
  }

  // Cont already contains the missing indices
  if(node_rank_){
    for(ULLInt in = 0; in < cont.size(); ++in){
      LLInt st_c = st[in];
      LLInt cont_c = cont[in];
      MatSetValues(ham_mat, 1, &cont_c, 1, &st_c, &ti, ADD_VALUES);
    }
  }

  MatAssemblyBegin(ham_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(ham_mat, MAT_FINAL_ASSEMBLY);

  MatSetOption(ham_mat, MAT_SYMMETRIC, PETSC_TRUE);
}
