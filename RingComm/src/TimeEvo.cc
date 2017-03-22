#include "TimeEvo.h"

TimeEvo::TimeEvo(const Environment &env)
{
  mpirank_ = env.mpirank;
  mpisize_ = env.mpisize;
  l_ = env.l;
  n_ = env.n;
}

TimeEvo::~TimeEvo()
{}

void TimeEvo::loschmidt_echo(const SparseOp &sparse, 
                             const unsigned int &iterations, 
                             const double *times, 
                             const double &tol, 
                             const int &maxits, 
                             const Vec &v, 
                             const Mat &ham_mat)
{
  MFN mfn;
  FN f;
  MFNCreate(PETSC_COMM_WORLD, &mfn);
  MFNSetOperator(mfn, ham_mat);
  MFNGetFN(mfn, &f);
  FNSetType(f, FNEXP);
  MFNSetTolerances(mfn, tol, maxits);

  MFNSetType(mfn, MFNEXPOKIT);
  MFNSetUp(mfn);

  MFNConvergedReason reason;

  PetscScalar l_echo;
  PetscReal loschmidt;

  Vec vec_help;
  VecDuplicate(v, &vec_help);
  VecCopy(v, vec_help);

  VecDot(v, v, &l_echo);
  loschmidt = (PetscRealPart(l_echo) * PetscRealPart(l_echo))
    + (PetscImaginaryPart(l_echo) * PetscImaginaryPart(l_echo));

  if(mpirank_ == 0){
      std::cout << "Time" << "\t" << "LE" << std::endl;
      std::cout << times[0] << "\t" << loschmidt << std::endl;
  }
  for(unsigned int tt = 1; tt < (iterations + 1); ++tt){

    FNSetScale(f, (times[tt] - times[tt - 1]) * PETSC_i, 1.0);
    MFNSolve(mfn, vec_help, vec_help);

    MFNGetConvergedReason(mfn, &reason);
    if(reason < 0){
      std::cerr << "Solver did not converge, aborting" << std::endl;
      std::cerr << "Change tolerance or maximum number of iterations" << std::endl;
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Echo
    VecDot(v, vec_help, &l_echo);
    loschmidt = (PetscRealPart(l_echo) * PetscRealPart(l_echo))
      + (PetscImaginaryPart(l_echo) * PetscImaginaryPart(l_echo));

    if(mpirank_ == 0){
        std::cout << times[tt] << "\t" << loschmidt << std::endl;
    }
  }
 
  MFNDestroy(&mfn);
  VecDestroy(&vec_help);
}


