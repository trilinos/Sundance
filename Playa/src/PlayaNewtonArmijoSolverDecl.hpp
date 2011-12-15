/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_NEWTON_ARMIJO_SOLVER_DECL_HPP
#define PLAYA_NEWTON_ARMIJO_SOLVER_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaNonlinearSolverBase.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * Playa implementation of Newton's method with Armijo line search.
 *
 * The solver's behavior is controlled by parameters in a ParameterList.
 * <ul>
 * <li> Scalar "Tau Relative" tolerance for relative error. Default value: 10 times machine epsilon
 * <li> Scalar "Tau Absolute" tolerance for absolute error. Default value 10 times machine epsilon
 * <li> Scalar "Alpha" constant in Armijo sufficient decrease condition. Default value 1.0e-4
 * <li> double "Step Reduction" factor by which to reduce step during line search. Default value: 0.5 
 * <li> int "Max Iterations" number of iterations to allow before failure. Default value 20.
 * <li> int "Max Backtracks" number of step reductions to allow before failure. Default value 20.
 * <li> int "Verbosity" amount of diagnostic output. Default value 0.
 * </ul>
 */
template <class Scalar>
class NewtonArmijoSolver : public NonlinearSolverBase<Scalar> 
{
public:
  /** */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  NewtonArmijoSolver(const ParameterList& params, 
    const LinearSolver<Scalar>& linSolver);

  /** */
  virtual ~NewtonArmijoSolver(){;}

  /** */
  SolverState<Scalar> solve(const NonlinearOperator<Scalar>& F,
    Vector<Scalar>& soln) const ;

  /* */
  GET_RCP(NonlinearSolverBase<Scalar>);

private:
  LinearSolver<Scalar> linSolver_;
  ScalarMag tauR_;
  ScalarMag tauA_;
  ScalarMag alpha_;
  ScalarMag stepReduction_;
  int maxIters_;
  int maxLineSearch_;
  int verb_;
    
};

  
}

#endif
