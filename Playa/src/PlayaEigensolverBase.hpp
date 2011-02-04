/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_EIGENSOLVERBASE_HPP
#define PLAYA_EIGENSOLVERBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp" 
#include "PlayaSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaLinearOperatorImpl.hpp"

namespace Playa
{
using Teuchos::ParameterList;

/**
 * Base class for eigensolvers for linear eigenvalue problems
 * \f[
 * K x = \lambda M x.
 * \f]
 */
template <class Scalar>
class EigensolverBase 
{
public:
  /** */
  EigensolverBase() : params_() {;}

  /** */
  EigensolverBase(const ParameterList& params) : params_(params) {;}

  /** */
  virtual ~EigensolverBase(){;}

  /** 
   * Solve a generalized eigensystem \f$K x = \lambda M x.\f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    const LinearOperator<Scalar>& M,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const = 0 ;

  /** 
   * Solve an eigensystem \f$K x = \lambda x.\f$
   */
  virtual void solve(
    const LinearOperator<Scalar>& K,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const 
    {
      LinearOperator<Scalar> M;
      solve(K,M,ev,ew);
    };

  /** 
   * Return the parameter list that was used to define this object. 
   */
  const ParameterList& params() const {return params_;}
  
private:
  ParameterList params_;
};  

}


#endif
