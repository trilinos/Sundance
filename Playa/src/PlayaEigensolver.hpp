/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_EIGENSOLVER_HPP
#define PLAYA_EIGENSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp" 
#include "PlayaSolverState.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaEigensolverBase.hpp"

namespace Playa
{
using Teuchos::ParameterList;

/**
 * Handle class for eigensolvers
 */
template <class Scalar>
class Eigensolver : public Playa::Handle<EigensolverBase<Scalar> >
{
public:
  /** */
  Eigensolver() : Handle<EigensolverBase<Scalar> >() {;}
  /** */
  Eigensolver( Playa::Handleable<EigensolverBase<Scalar> >* rawPtr) 
    : Handle<EigensolverBase<Scalar> >(rawPtr) {;}
  /** */
  Eigensolver(const RCP<EigensolverBase<Scalar> >& smartPtr)
    : Handle<EigensolverBase<Scalar> >(smartPtr) {;}
  

  /** */
  void solve(
    const LinearOperator<Scalar>& K,
    const LinearOperator<Scalar>& M,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const 
    {
      this->ptr()->solve(K, M, ev, ew);
    }

  /** */
  void solve(
    const LinearOperator<Scalar>& K,
    Array<Vector<Scalar> >& ev,
    Array<std::complex<Scalar> >& ew) const 
    {
      this->ptr()->solve(K, ev, ew);
    }

  /** */
  const ParameterList& params() const 
    {
      TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Eigensolver::parameters()");
      return this->ptr()->params();
    }
  
  /** */
  ParameterList& params() 
    {
      TEST_FOR_EXCEPTION(this->ptr().get()==0, std::runtime_error,
                       "null pointer in Eigensolver::parameters()");
      return this->ptr()->params();
    }
  

  static FancyOStream& os()
    {
      return Out::os();
    }
};


  
}


#endif
