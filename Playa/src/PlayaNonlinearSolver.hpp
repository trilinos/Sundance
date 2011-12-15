/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_NONLINEAR_SOLVER_HPP
#define PLAYA_NONLINEAR_SOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaNonlinearSolverBase.hpp"
#include "PlayaHandle.hpp"

namespace Playa
{

/**
 *
 */
template <class Scalar>
class NonlinearSolver : public Handle<NonlinearSolverBase<Scalar> >
{
public:
  /* boilerplate ctors */
  HANDLE_CTORS(NonlinearSolver<Scalar>, NonlinearSolverBase<Scalar>);

  /** */
  SolverState<Scalar> solve(const NonlinearOperator<Scalar>& F,
    Vector<Scalar>& soln) const {return this->ptr()->solve(F, soln);}
};

  
}

#endif
