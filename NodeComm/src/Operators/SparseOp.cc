#include "SparseOp.h"

/*******************************************************************************/
// Single custom constructor for this class.
// Creates the Hamiltonian matrix depending on the basis chosen.
/*******************************************************************************/
SparseOpNC::SparseOpNC(const EnvironmentNC &env, const BasisNC &basis)
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

  MatCreate(PETSC_COMM_WORLD, &HamMat);
  MatSetSizes(HamMat, nlocal_, nlocal_, basis_size_, basis_size_);
  MatSetType(HamMat, MATMPIAIJ);
}

/*******************************************************************************/
// Copy constructor
/*******************************************************************************/
SparseOpNC::SparseOpNC(const SparseOpNC &rhs)
{
  std::cout << "Copy constructor (ham matrix) has been called!" << std::endl;

  l_ = rhs.l_;
  n_ = rhs.n_;
  mpirank_ = rhs.mpirank_;
  mpisize_ = rhs.mpisize_;
  node_rank_ = rhs.node_rank_;
  node_size_ = rhs.node_size_;
  nlocal_ = rhs.nlocal_;
  start_ = rhs.start_;
  end_ = rhs.end_;
  basis_size_ = rhs.basis_size_;
  
  MPI_Comm_dup(rhs.node_comm_, &node_comm_);
  MatDuplicate(rhs.HamMat, MAT_COPY_VALUES, &HamMat);
}

/*******************************************************************************/
// Assignment operator
/*******************************************************************************/
SparseOpNC &SparseOpNC::operator=(const SparseOpNC &rhs)
{
  std::cout << "Assignment operator (ham matrix) has been called!" << std::endl;
    
  if(this != &rhs){
    MatDestroy(&HamMat);
    MPI_Comm_free(&node_comm_);    

    l_ = rhs.l_;
    n_ = rhs.n_;
    mpirank_ = rhs.mpirank_;
    mpisize_ = rhs.mpisize_;
    node_rank_ = rhs.node_rank_;
    node_size_ = rhs.node_size_;
    nlocal_ = rhs.nlocal_;
    start_ = rhs.start_;
    end_ = rhs.end_;
    basis_size_ = rhs.basis_size_;
  
    MPI_Comm_dup(rhs.node_comm_, &node_comm_);
    MatDuplicate(rhs.HamMat, MAT_COPY_VALUES, &HamMat);
  }

  return *this;
}

SparseOpNC::~SparseOpNC()
{
  MatDestroy(&HamMat);
  MPI_Comm_free(&node_comm_);
}

/*******************************************************************************/
// Determines the sparsity pattern to allocate memory only for the non-zero 
// entries of the matrix
/*******************************************************************************/
void SparseOpNC::determine_allocation_details_(LLInt *int_basis, 
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

          LLInt new_int1 = UtilsNC::binary_to_int(bitset, l_);
          // Loop over all states and look for a match
          LLInt match_ind1;
          if(node_rank_){
            match_ind1 = UtilsNC::binsearch(int_basis, nlocal_, new_int1); 
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
            match_ind1 = UtilsNC::binsearch(int_basis, basis_size_, new_int1);
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

          LLInt new_int0 = UtilsNC::binary_to_int(bitset, l_);
          // Loop over all states and look for a match
          LLInt match_ind0;
          if(node_rank_){
            match_ind0 = UtilsNC::binsearch(int_basis, nlocal_, new_int0); 
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
            match_ind0 = UtilsNC::binsearch(int_basis, basis_size_, new_int0);
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
        LLInt m_ind = UtilsNC::binsearch(int_basis, basis_size_, cont[ii]);
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
void SparseOpNC::construct_AA_hamiltonian(LLInt *int_basis, 
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
  MatMPIAIJSetPreallocation(HamMat, 0, d_nnz, 0, o_nnz);

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
        MatSetValues(HamMat, 1, &state, 1, &state, &osc_term, ADD_VALUES);
        
        int next_site1 = (site + 1) % l_;

        // If there's a particle in next site, do nothing
        if(bitset[next_site1] == 1){
          // Accumulate 'V' terms
          MatSetValues(HamMat, 1, &state, 1, &state, &Vi, ADD_VALUES);
          continue;
        }
        // Otherwise do a swap
        else{
          bitset[next_site1] = 1;
          bitset[site]       = 0;

          LLInt new_int1 = UtilsNC::binary_to_int(bitset, l_);
          // Loop over all states and look for a match
          LLInt match_ind1;
          if(node_rank_){
            match_ind1 = UtilsNC::binsearch(int_basis, nlocal_, new_int1); 
            if(match_ind1 == -1){
              continue;
            }
            else{
              match_ind1 += start_;
              MatSetValues(HamMat, 1, &match_ind1, 1, &state, &ti, ADD_VALUES);
            }
          }
          else{
            match_ind1 = UtilsNC::binsearch(int_basis, basis_size_, new_int1);
            MatSetValues(HamMat, 1, &match_ind1, 1, &state, &ti, ADD_VALUES);
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

          LLInt new_int0 = UtilsNC::binary_to_int(bitset, l_);
          // Loop over all states and look for a match
          LLInt match_ind0;
          if(node_rank_){
            match_ind0 = UtilsNC::binsearch(int_basis, nlocal_, new_int0); 
            if(match_ind0 == -1){
              continue;
            }
            else{
              match_ind0 += start_;
              MatSetValues(HamMat, 1, &match_ind0, 1, &state, &ti, ADD_VALUES);
            }
          }
          else{
            match_ind0 = UtilsNC::binsearch(int_basis, basis_size_, new_int0);
            MatSetValues(HamMat, 1, &match_ind0, 1, &state, &ti, ADD_VALUES);
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
      MatSetValues(HamMat, 1, &cont_c, 1, &st_c, &ti, ADD_VALUES);
    }
  }

  MatAssemblyBegin(HamMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(HamMat, MAT_FINAL_ASSEMBLY);

  MatSetOption(HamMat, MAT_SYMMETRIC, PETSC_TRUE);
}
