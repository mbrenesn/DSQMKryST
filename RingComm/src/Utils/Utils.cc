#include "Utils.h"

namespace Utils
{
  LLInt mod(LLInt a, LLInt b)
  {
    return (a % b + b) % b;
  }

  LLInt binary_to_int(boost::dynamic_bitset<> bs, unsigned int l)
  {
    LLInt integer = 0;

    for(unsigned int i = 0; i < l; ++i){
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
  LLInt binsearch(const LLInt *array, LLInt len, LLInt value)
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
  LLInt get_neel_index(const Environment &env, const Basis &bas)
  {
    LLInt index;
    if(env.l / 2 != env.n){
      std::cerr << "Not implemented!" << std::endl;
      std::cerr << "Neel state has only been implemented for half-filled systems" << std::endl;
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
  
    boost::dynamic_bitset<> neel(env.l, 1);
    for(unsigned int site = 0; site < env.l; site += 2){
      neel.set(site);
    }
  
    LLInt neel_int = binary_to_int(neel, env.l);
    index = binsearch(bas.int_basis, bas.basis_local, neel_int);
    if(index != -1) index += bas.start;
  
    return index;
  }
  
  /*******************************************************************************/
  // Pick initial index out of the basis randomly 
  /*******************************************************************************/
  LLInt get_random_index(const Environment &env, const Basis &bas, bool wtime, bool verbose)
  {
    LLInt pick_ind;
    boost::random::mt19937 gen;
  
    if(wtime) gen.seed(static_cast<LLInt>(std::time(0)));
  
    if(env.mpirank == 0){
      boost::random::uniform_int_distribution<LLInt> dist(0, bas.basis_size - 1);
      pick_ind = dist(gen);
    }
    MPI_Bcast(&pick_ind, 1, MPI_LONG_LONG_INT, 0, PETSC_COMM_WORLD);
  
    bool check = false;
    if(pick_ind >= bas.start && pick_ind < bas.end) check = true;
      
    if(verbose && check){
      std::cout << "Initial state randomly chosen: " << bas.int_basis[pick_ind] << std::endl;
      std::cout << "With binary representation: " << std::endl;
      boost::dynamic_bitset<> bs(env.l, bas.int_basis[pick_ind - bas.start]);
      std::cout << bs << std::endl;
    }
  
    return pick_ind;
  }
}
