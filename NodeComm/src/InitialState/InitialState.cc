#include "InitialState.h"

/*******************************************************************************/
// Single custom constructor for this class.
// Creates the initial state object.
/*******************************************************************************/
InitialStateNC::InitialStateNC(const EnvironmentNC &env,
                               const BasisNC &basis)
{
  l_ = env.l;
  n_ = env.n;
  mpirank_ = env.mpirank;
  mpisize_ = env.mpisize;
  nlocal_ = basis.nlocal;
  start_ = basis.start;
  end_ = basis.end;
  basis_size_ = basis.basis_size;

  VecCreateMPI(PETSC_COMM_WORLD, nlocal_, basis_size_, &InitialVec);
}

/*******************************************************************************/
// Copy constructor
/*******************************************************************************/
InitialStateNC::InitialStateNC(const InitialStateNC &rhs)
{
  std::cout << "Copy constructor (initial state) has been called!" << std::endl;

  l_ = rhs.l_;
  n_ = rhs.n_;
  mpirank_ = rhs.mpirank_;
  mpisize_ = rhs.mpisize_;
  nlocal_ = rhs.nlocal_;
  start_ = rhs.start_;
  end_ = rhs.end_;
  basis_size_ = rhs.basis_size_;

  VecDuplicate(rhs.InitialVec, &InitialVec);
  VecCopy(rhs.InitialVec, InitialVec);
}

/*******************************************************************************/
// Assignment operator
/*******************************************************************************/
InitialStateNC &InitialStateNC::operator=(const InitialStateNC &rhs)
{
  std::cout << "Assignment operator (diagonal op) has been called!" << std::endl;
    
  if(this != &rhs){
    VecDestroy(&InitialVec);
      
    l_ = rhs.l_;
    n_ = rhs.n_;
    mpirank_ = rhs.mpirank_;
    mpisize_ = rhs.mpisize_;
    nlocal_ = rhs.nlocal_;
    start_ = rhs.start_;
    end_ = rhs.end_;
    basis_size_ = rhs.basis_size_;

    VecDuplicate(rhs.InitialVec, &InitialVec);
    VecCopy(rhs.InitialVec, InitialVec);
  }

  return *this;
}

InitialStateNC::~InitialStateNC()
{
  VecDestroy(&InitialVec);
}

/*******************************************************************************/
// Neel state
/*******************************************************************************/
void InitialStateNC::neel_initial_state(LLInt *int_basis)
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

  if(mpirank_ == 0){
    LLInt neel_int = UtilsNC::binary_to_int(neel, l_);
    index = UtilsNC::binsearch(int_basis, nlocal_, neel_int);
    VecSetValue(InitialVec, index, 1.0, INSERT_VALUES);
  }

  VecAssemblyBegin(InitialVec);
  VecAssemblyEnd(InitialVec);
}

/*******************************************************************************/
// Initial random state out of the computational basis
/*******************************************************************************/
void InitialStateNC::random_initial_state(LLInt *int_basis,
                                          bool wtime,
                                          bool verbose)
{
  LLInt pick_ind;
  boost::random::mt19937 gen;

  if(wtime) gen.seed(static_cast<LLInt>(std::time(0)));

  if(mpirank_ == 0){
    boost::random::uniform_int_distribution<LLInt> dist(0, basis_size_ - 1);
    pick_ind = dist(gen);

    if(verbose){
      std::cout << "Initial state randomly chosen: " << int_basis[pick_ind] << std::endl;
      std::cout << "With binary representation: " << std::endl;
      boost::dynamic_bitset<> bs(l_, int_basis[pick_ind - start_]);
      std::cout << bs << std::endl;
    }
    VecSetValue(InitialVec, pick_ind, 1.0, INSERT_VALUES);
  }

  VecAssemblyBegin(InitialVec);
  VecAssemblyEnd(InitialVec);
}
