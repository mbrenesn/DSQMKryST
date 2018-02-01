/** @addtogroup NodeComm
 * @{
 */
/**
 * \class EnvironmentNC.
 * \ingroup NodeComm
 * \brief Initializes the MPI, PETSc and SLEPc environments.
 *
 * This class should be instantiated at the beginning of the program and at the end
 * of execution the environment will close automatically
 * Refer to Section 3 of the manuscript in /docs.
 * For the particular case of the Node communication approach, the basis elements are 
 * fully distributed among all available processing elements, with the exception of the first process
 * of each computational node; which allocates, computes and holds all the elements of the basis. 
 */
#ifndef __ENVIRONMENT_H
#define __ENVIRONMENT_H

#include <iostream>
#include <stdint.h>

#include <petscsys.h>
#include <slepcmfn.h>
#include <slepceps.h>

typedef unsigned long long ULLInt;
typedef PetscInt LLInt;

class EnvironmentNC
{
  public:
    /** \brief Creates an instance of class Environment.
      * \param argc Required to parse PETSc/SLEPc's options.
      * \param argv Required to parse PETSc/SLEPc's options.
      * \param l Number of sites.
      * \param n Subspace descriptor (number of particles).
      *
      * This is the only available constructor of this class. This constructor is used to 
      * initialise PETSc, SLEPc and MPI environments and should be instantiated at the 
      * beginning of the program.
      */
    EnvironmentNC(int argc, 
                char **argv, 
                unsigned int l, 
                unsigned int n);
    /** \brief Destructor.
      * 
      * Closes the MPI, PETSc and SLEPc environment. Destroys the node communicator.
      */
    ~EnvironmentNC();
    /** \brief Computes the dimension of the Hilbert space.
      */
    LLInt basis_size() const;
    /** \brief Computes the section and global indices of locally owned elements.
      * \param b_size Dimension of the Hilbert space.
      * \param nlocal Local section of the process (amount of elements/number of rows).
      * \param start Global index refering to the local process.
      * \param end Global index refering to the local process.
      * 
      * This parallel distribution agrees with PETSc's row-wise distribution
      * Refer to Algorithm 1 of the manuscript in /docs.
      */
    void distribution(PetscInt b_size, 
                      PetscInt &nlocal, 
                      PetscInt &start, 
                      PetscInt &end) const;
    unsigned int l; ///< Number of sites.
    unsigned int n; ///< Subspace descriptor (number of particles).
    PetscMPIInt mpirank; ///< Index of the local processor.
    PetscMPIInt mpisize; ///< Total number of processors.
    PetscMPIInt node_rank; ///< Rank respective to the node
    PetscMPIInt node_size; ///< Number of processes per node
    MPI_Comm node_comm; ///< The MPI communicator respective of the node
  private:
};
#endif
/** @}*/
