#include "Basis.h"

/*******************************************************************************/
// Custom/only constructor
/*******************************************************************************/
Basis::Basis(const Environment &env)
{
  l_ = env.l;
  n_ = env.n;
  basis_size = env.basis_size();
  env.distribution(basis_size, nlocal, start, end);

  basis_local = nlocal;
  basis_start = start;

  int_basis = new LLInt[basis_local];
}

/*******************************************************************************/
// Copy constructor
/*******************************************************************************/
Basis::Basis(const Basis &rhs)
{
  std::cout << "Copy constructor (basis) has been called!" << std::endl;

  l_ = rhs.l_;
  n_ = rhs.n_;
  basis_size = rhs.basis_size;
  nlocal = rhs.nlocal;
  start = rhs.start;
  end = rhs.end;
  basis_local = rhs.basis_local;
  basis_start = rhs.basis_start;

  int_basis = new LLInt[basis_local];
  for(LLInt i = 0; i < basis_local; ++i)
    int_basis[i] = rhs.int_basis[i];
}

/*******************************************************************************/
// Assignment operator
/*******************************************************************************/
Basis &Basis::operator=(const Basis &rhs)
{
  std::cout << "Assignment operator (basis) has been called!" << std::endl;

  delete [] int_basis;
  l_ = rhs.l_;
  n_ = rhs.n_;
  basis_size = rhs.basis_size;
  nlocal = rhs.nlocal;
  start = rhs.start;
  end = rhs.end;
  basis_local = rhs.basis_local;
  basis_start = rhs.basis_start;

  int_basis = new LLInt[basis_local];
  for(LLInt i = 0; i < basis_local; ++i)
    int_basis[i] = rhs.int_basis[i];

  return *this;
}

Basis::~Basis()
{
  delete [] int_basis;
}

/*******************************************************************************/
// Helper function.
/*******************************************************************************/
LLInt Basis::factorial_(LLInt n)
{
  return (n == 1 || n == 0) ? 1 : factorial_(n - 1) * n;
}

/*******************************************************************************/
// Returns the smallest possible integer that can be expressed with a given
// binary combination.
/*******************************************************************************/
LLInt Basis::first_int_()
{
  LLInt first = 0;
  for(LLInt i = 0; i < n_; ++i){
    first += 1 << i;
  }

  LLInt w = first;
  for(LLInt i = 0; i < basis_start; ++i){
    LLInt t = (first | (first - 1)) + 1;
    w = t | ((((t & -t) / (first & -first)) >> 1) - 1);
    first = w;
  }

  return first;
}

/*******************************************************************************/
// Computes next bit permutation lexicographically, so for a given value of
// smallest ineteger, computes the next lowest integer that has a different 
// bit combination. Since we know the number of combinations possible, this 
// returns an array which contains all possible combinations represented as
// integer values.
/*******************************************************************************/
void Basis::construct_int_basis()
{
  LLInt w;                                         // Next permutation of bits
  LLInt first = first_int_();

  int_basis[0] = first;

  for(LLInt i = 1; i < basis_local; ++i){
    LLInt t = (first | (first - 1)) + 1;
    w = t | ((((t & -t) / (first & -first)) >> 1) - 1);
    
    int_basis[i] = w;

    first = w;
  }
}

/*******************************************************************************/
// Print to std out
/*******************************************************************************/
void Basis::print_basis(const Environment &env, 
                        bool bits)
{
  std::cout << "Global rank: " << env.mpirank << std::endl;
  for(LLInt i = 0; i < basis_local; ++i){
    if(bits){
      boost::dynamic_bitset<> bs(l_, int_basis[i]);
      std::cout << bs << std::endl;
    }
    else
      std::cout << int_basis[i] << std::endl;
  }
}

/*******************************************************************************/
// Given an integer represented basis we now compute the basis in binary no
// tation, boost library provides best results on the long run.
// See Pieterse, et al. (2010) for a reference
/*******************************************************************************/
void Basis::construct_bit_basis(boost::dynamic_bitset<> *bit_basis)
{
  for(LLInt i = 0; i < basis_local; ++i){
    boost::dynamic_bitset<> bs(l_, int_basis[i]);
    bit_basis[i] = bs;
  }
}
