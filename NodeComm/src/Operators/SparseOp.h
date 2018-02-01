/** @addtogroup NodeComm
 * @{
 */
/**
 * \class SparseOpNC.
 * \ingroup NodeComm
 * \brief Related to the matrix representation of the Hamiltonian of the quantum system.
 *
 * The Hamiltonian matrix itself is a public member of this class and is row-wise distributed among
 * processing elements. For the Node communicator approach, the row-wise distribution of the matrix
 * is retained, the only difference is the way the matrix elements are computed. See Section 3.1 
 * (Node communicator) of the manuscript in /docs for more details.
 */
#ifndef __SPARSEOP_H
#define __SPARSEOP_H

#include "../Environment/Environment.h"
#include "../Utils/Utils.h"
#include "../Basis/Basis.h"

class SparseOpNC
{
  public:
    /** \brief Creates an instance of class SparseOp.
      * \param env An instance of the class Environment.
      * \param basis An instance of the class Basis.
      *
      * This is the only available constructor of this class. After creating an instance of this class
      * one can call the construct_AA_hamiltonian(...) method to introduce parameters into the matrix
      * The matrix itself is a public member of the class, row-wise distributed across processing 
      * elements. Please refer to Figure 1 and Section 3 (Hamiltonian matrix construction) of the 
      * manuscript in /docs for more information.
      */
    SparseOpNC(const EnvironmentNC &env, 
               const BasisNC &basis);
    /** \brief Destructor.
      * 
      * Deallocates and destroys the Hamiltonian matrix, no need to call MatDestroy() on the matrix. 
      */ 
    ~SparseOpNC();
    /// Copy constructor.
    SparseOpNC(const SparseOpNC &rhs);
    /// Overloading of the assignment operator.
    SparseOpNC &operator=(const SparseOpNC &rhs);
    /** \brief Allocates memory and inserts elements to Hamiltonian matrix.
      * 
      * This should be called after creating an instance of SparseOp and before using time-evolution
      * routines. The member HamMat is a matrix of type MATMPIAIJ from PETSc, for which memory is 
      * preallocated, distributed and elements are added by this routine. The main communication
      * described in Algorithm 5 and Section 3.1 (node communicator approach) in the manuscript
      * located in /docs is used in this routine.
      */
    void construct_AA_hamiltonian(LLInt *int_basis, 
                                  double V,
                                  double t, 
                                  double h,
                                  double beta);
    Mat HamMat; ///< The Hamiltonian matrix, row-wise distributed. PETSc MATMPIAIJ object.

  private:
    unsigned int l_; ///< Number of sites.
    unsigned int n_; ///< Subspace descriptor.
    PetscMPIInt mpirank_; ///< Index of the local processor.
    PetscMPIInt mpisize_; ///< Total number of processors.
    PetscMPIInt node_rank_; ///< Local rank to specific node
    PetscMPIInt node_size_; ///< Number of processing elements in node
    MPI_Comm node_comm_; ///< The node MPI communicator
    LLInt basis_size_; ///< Dimension of the Hilbert space.
    PetscInt nlocal_; ///< Local amount of rows owned by processor (PETSc).
    PetscInt start_; ///< Global index (PETSc).
    PetscInt end_; ///< Global index (PETSc).
    /** \brief A communication routine, wrapper to MPI_Allgather.
      * 
      * Section 3.1 and Algorithm 5 of the document in /docs for more details 
      */
    void gather_nonlocal_values_(LLInt *start_inds);
    /** \brief Computes the number of non-zero elements.
      * 
      * For good performance, the sparse matrix that represents the Hamiltonian of the system
      * has to be preallocated in memory. This routine is called internally by contruct_AA_hamiltonian()
      * to allocate memory for the matrix. The main communication procedure described in Section 3.1 
      * (node communicator approach) and Algorithm 5 of the manuscript in /docs is implemented here. 
      */
    void determine_allocation_details_(LLInt *int_basis, 
                                       std::vector<LLInt> &cont,
                                       std::vector<LLInt> &st, 
                                       PetscInt *diag, 
                                       PetscInt *off);
};
#endif
/** @}*/
