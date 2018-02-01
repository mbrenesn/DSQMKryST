#include "InitialState.h"

/*******************************************************************************/
// Single custom constructor for this class.
// Creates the initial state object.
/*******************************************************************************/
InitialStateRC::InitialStateRC(const EnvironmentRC &env,
                               const BasisRC &basis)
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
InitialStateRC::InitialStateRC(const InitialStateRC &rhs)
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
InitialStateRC &InitialStateRC::operator=(const InitialStateRC &rhs)
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

InitialStateRC::~InitialStateRC()
{
  VecDestroy(&InitialVec);
}

/*******************************************************************************/
// Neel state
/*******************************************************************************/
void InitialStateRC::neel_initial_state(LLInt *int_basis)
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

  LLInt neel_int = UtilsRC::binary_to_int(neel, l_);
  index = UtilsRC::binsearch(int_basis, nlocal_, neel_int);
  if(index != -1){
    index += start_;
    VecSetValue(InitialVec, index, 1.0, INSERT_VALUES);
  }
  VecAssemblyBegin(InitialVec);
  VecAssemblyEnd(InitialVec);
}

/*******************************************************************************/
// Initial random state out of the computational basis
/*******************************************************************************/
void InitialStateRC::random_initial_state(LLInt *int_basis,
                                          bool wtime,
                                          bool verbose)
{
  LLInt pick_ind;
  boost::random::mt19937 gen;

  if(wtime) gen.seed(static_cast<LLInt>(std::time(0)));

  if(mpirank_ == 0){
    boost::random::uniform_int_distribution<LLInt> dist(0, basis_size_ - 1);
    pick_ind = dist(gen);
  }
  MPI_Bcast(&pick_ind, 1, MPI_LONG_LONG_INT, 0, PETSC_COMM_WORLD);

  bool check = false;
  if(pick_ind >= start_ && pick_ind < end_) check = true;
    
  if(verbose && check){
    std::cout << "Initial state randomly chosen: " << int_basis[pick_ind] << std::endl;
    std::cout << "With binary representation: " << std::endl;
    boost::dynamic_bitset<> bs(l_, int_basis[pick_ind - start_]);
    std::cout << bs << std::endl;
  }
  VecSetValue(InitialVec, pick_ind, 1.0, INSERT_VALUES);
  VecAssemblyBegin(InitialVec);
  VecAssemblyEnd(InitialVec);
}
