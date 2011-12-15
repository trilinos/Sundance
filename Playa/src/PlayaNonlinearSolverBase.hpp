/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_NONLINEARSOLVERBASE_HPP
#define PLAYA_NONLINEARSOLVERBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"

namespace Playa
{
using namespace Teuchos;
template<class Scalar>
class NonlinearOperator;

/**
 * Base interface for nonlinear solvers
 */
template <class Scalar>
class NonlinearSolverBase : public Handleable<NonlinearSolverBase<Scalar> >
{
public:
  /** */
  NonlinearSolverBase(const ParameterList& params = ParameterList());

  /** */
  virtual ~NonlinearSolverBase(){;}

  /** */
  virtual SolverState<Scalar> solve(const NonlinearOperator<Scalar>& F,
    Vector<Scalar>& soln) const = 0  ;

protected:

  const ParameterList& params() const {return params_;}

private:
  ParameterList params_;
};

  
template <class Scalar> inline
NonlinearSolverBase<Scalar>
::NonlinearSolverBase(const ParameterList& params)
  : params_(params)
{;}
  
}

#endif
