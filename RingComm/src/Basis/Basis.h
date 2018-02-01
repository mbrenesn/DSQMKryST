/** @addtogroup RingComm
 * @{
 */
/**
 * \class BasisRC
 * \ingroup RingComm
 * \brief A computational representation of the Hilbert space basis.
 *        Refer to Section 3 of the manuscript in /docs.
 *
 * For the particular case of the Ring communication approach, the basis integer representations are 
 * fully distributed among all available processing elements.
 */
#ifndef __BASIS_H
#define __BASIS_H

#include <boost/dynamic_bitset.hpp>

#include "../Environment/Environment.h"

class BasisRC
{
  public:
    /** \brief Creates an instance of class Basis.
      * \param env An instance of the class Environment.
      *
      * This is the only available constructor of this class. A basis is required to provide a matrix
      * representation of a Hamiltonian operator and carry out other operations in quantum mechanics. 
      * We use an integer representation of the states in the basis of the Hilbert space of the system, 
      * please refer to Section 2 (Background) and Section 3 (Basis representation) of the manuscript
      * in /docs.
      */
    BasisRC(const EnvironmentRC &env);
    /** \brief Destructor.
      * 
      * Deallocates memory needed for the computational representation of the basis.
      */ 
    ~BasisRC();
    /// Copy constructor.
    BasisRC(const BasisRC &rhs);
    /// Overloading of the assignment operator.
    BasisRC &operator=(const BasisRC &rhs);
    /** \brief Computes and distributes (using MPI) the representation of the basis elements.
      *
      * This should be called after creating an instance of Basis and before
      * constructing a Hamiltonian matrix. For more details, refer to Section 3.1 (Ring exchange
      * approach) of the manuscript in /docs.
      */
    void construct_int_basis();
    /** \brief Outputs the integer representation of the basis to stdout.
      * \param env An instance of class Environment.
      * \param bits If bits = true, outputs a bit representation.
      */
    void print_basis(const EnvironmentRC &env, 
                     bool bits = false);
    /** \brief Computes and distributes (using MPI) the representation of the basis elements
      *        in binary form, this is normally used for visualisation purposes.
      */
    void construct_bit_basis(boost::dynamic_bitset<> *bit_basis);
    LLInt basis_size; ///< Dimension of the Hilbert space.
    PetscInt basis_local; ///< Local value (MPI) of the dimension of the Hilbert space.
    PetscInt basis_start; ///< Global index per processor.
    PetscInt nlocal; ///< Local amount of rows owned by processor (PETSc).
    PetscInt start; ///< Global index (PETSc).
    PetscInt end; ///< Global index (PETSc).
    LLInt *int_basis; ///< Container of the elements of the basis, locally owned.
  
  private:
    unsigned int l_; ///< Number of sites.
    unsigned int n_; ///< Subspace descriptor.
    /** \brief Computes the factorial of an integer.
     *  \param n Integer value.
     *  \return The factorial of the number.
     */
    LLInt factorial_(LLInt n); 
    /** \brief Computes the first integer representation of the basis vector locally owned.
     *  \return The integer representation.
     */
    LLInt first_int_();
};
#endif
/** @}*/
