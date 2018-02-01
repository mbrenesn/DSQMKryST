/** @addtogroup RingComm
 * @{
 */
/**
 * \class KrylovEvoRC.
 * \ingroup RingComm
 * \brief This class is used to evaluate the dynamics of the system using the Krylov approach.
 *
 * The method used to evaluate the action of the propagator is described in Section 2 and Section 3
 * (Time evolution) of the manuscript in /docs. The Krylov subspace method applied in this approach
 * is the Arnoldi method in conjuction with irreducible Pade approximations (Figure 1 of the manuscript).
 * The method used in this implementation is the one implemented in SLEPc as part of the MFN section,
 * which is the same method described by R. Sidje (https://www.maths.uq.edu.au/expokit/)
 * For more information refer to the manuscript in /docs and the SLEPc manual: 
 * http://slepc.upv.es/documentation/
 */
#ifndef __KRYLOV_EVO_H
#define __KRYLOV_EVO_H

#include "../Environment/Environment.h"
#include "../Basis/Basis.h"

class KrylovEvoRC
{
  public:
    /** \brief Creates an instance of class KrylovEvo.
      * \param ham_mat The Hamiltonian matrix, or any other matrix object from PETSc.
      * \param tol Tolerance of the algorithm.
      * \param max_kryt_its Maximum amount of iterations of the algorithm.
      *
      * This is the only available constructor of this class. After creating an instance of this class
      * the MFN and FN environments are set with the parameters given.
      */
    KrylovEvoRC(const Mat &ham_mat,
                const double &tol,
                const int &max_kryt_its);
    /** \brief Destructor.
      * 
      * Deallocates and destroys objects associated with the MFN component of SLEPc.
      */ 
    ~KrylovEvoRC();
    MFNConvergedReason reason; ///< Object related to the convergence of the algorithm.
                               ///< If !=0, then the algorithm failed to converge with given parameters.
    /** \brief Time evolution routine.
      * \param final_time Final time value.
      * \param initial_time Initial time value.
      * \param vec A vector that represents the initial state at time = initial_time.
      *
      * This method is used to evolve in time a given state. The state is replaced with it's time-evolved
      * counterpart.
      */ 
    void krylov_evo(const double &final_time,
                    const double &initial_time,
                    Vec &vec);
  
  private:
    MFN mfn_; ///< MFN component object, containing details related to parameters of the algorithm.
    FN f_; ///< FN component object, containing details related to the function to be applied to the
           ///< operator, exponential in this particular case.
};
#endif
/** @}*/
