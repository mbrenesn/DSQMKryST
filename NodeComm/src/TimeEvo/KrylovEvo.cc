#include "KrylovEvo.h"

KrylovEvoNC::KrylovEvoNC(const Mat &ham_mat,
                         const double &tol,
                         const int &max_kryt_its)
{
  MFNCreate(PETSC_COMM_WORLD, &mfn_);
  MFNSetOperator(mfn_, ham_mat);
  MFNGetFN(mfn_, &f_);
  FNSetType(f_, FNEXP);
  MFNSetTolerances(mfn_, tol, max_kryt_its);

  MFNSetType(mfn_, MFNEXPOKIT);
  MFNSetUp(mfn_);
}

KrylovEvoNC::~KrylovEvoNC()
{
  MFNDestroy(&mfn_);
}

void KrylovEvoNC::krylov_evo(const double &final_time,
                             const double &initial_time,
                             Vec &vec)
{
  FNSetScale(f_, (final_time - initial_time) * PETSC_i, 1.0);
  MFNSolve(mfn_, vec, vec);

  MFNGetConvergedReason(mfn_, &reason);

  if(reason < 0){
    std::cerr << "Krylov solver did not converge, aborting" << std::endl;
    std::cerr << "Change tolerance or maximum number of iterations" << std::endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
}
