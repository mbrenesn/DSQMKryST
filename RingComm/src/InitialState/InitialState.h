/**
 * \class InitialState.
 *
 * \brief A class to construct initial states for time evolution.
 *
 * Currently only supports two different initial states: a Neel state and a random state
 * We note that a random state in this context refers to an initial state randomly picked
 * out of the computational basis.
 */
#ifndef __INITIAL_STATE_H
#define __INITIAL_STATE_H

#include "../Environment/Environment.h"
#include "../Utils/Utils.h"
#include "../Basis/Basis.h"

class InitialState
{
  public:
    /** \brief Creates an instance of class InitialState.
      * \param env An instance of the class Environment.
      * \param basis An instance of class Basis.
      *
      * This is the only available constructor of this class.
      */
    InitialState(const Environment &env,
                 const Basis &basis);
    /** \brief Destructor.
      * 
      * Destroys the initial state vector automatically.
      */ 
    ~InitialState();
    /// Copy constructor.
    InitialState(const InitialState &rhs);
    /// Overloading of the assignment operator.
    InitialState &operator=(const InitialState &rhs);
    Vec InitialVec; ///< Initial state represented as a vector in Hilbert space with same parallel layout.
    /** \brief Method to compute the Neel state.
      * \param int_basis The integer basis, a member of class Basis.
      */     
    void neel_initial_state(LLInt *int_basis);
    /** \brief Method to compute a random initial state.
      * \param int_basis The integer basis, a member of class Basis.
      * \param wtime If true, random state changes with each execution based on current time.
      * \param verbose If true, prints to stdout the random state chosen.
      *
      * RNG is Mersenne-Twister from Boost.
      */     
    void random_initial_state(LLInt *int_basis,
                              bool wtime = false,
                              bool verbose = false);
  private:
    unsigned int l_; ///< Number of sites.  
    unsigned int n_; ///< Subspace descriptor.
    PetscMPIInt mpirank_; ///< Index of the local processor.
    PetscMPIInt mpisize_; ///< Total number of processors.
    LLInt basis_size_; ///< Dimension of the Hilbert space. 
    PetscInt nlocal_; ///< Local amount of rows owned by processor (PETSc).
    PetscInt start_; ///< Global index (PETSc).
    PetscInt end_; ///< Global index (PETSc).
};
#endif
