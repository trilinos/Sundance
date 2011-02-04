/* @HEADER@ */
/* @HEADER@ */

#ifndef PLAYA_NONLINEAROPERATOR_HPP
#define PLAYA_NONLINEAROPERATOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaNonlinearOperatorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Playa
{
using namespace Teuchos;

/** 
 * User-level nonlinear operator class
 */
template <class Scalar>
class NonlinearOperator : public Handle<NonlinearOperatorBase<Scalar> >
{
public:
  /* boilerplate ctors */
  HANDLE_CTORS(NonlinearOperator<Scalar>, NonlinearOperatorBase<Scalar>);

  /** */
  VectorSpace<Scalar> domain() const 
    {return this->ptr()->domain();}

  /** */
  VectorSpace<Scalar>  range() const 
    {return this->ptr()->range();}

  /** */
  void setEvalPt(const Vector<double>& evalPt)
    {
      this->ptr()->setEvalPt(evalPt);
    }
      
  /** */
  LinearOperator<Scalar> getJacobian() const 
    {
      return this->ptr()->getJacobian();
    }

  /** */
  Vector<double> getFunctionValue() const 
    {
      return this->ptr()->getFunctionValue();
    }

      

  /** */
  Vector<double> getInitialGuess() const 
    {
      return this->ptr()->getInitialGuess();
    }

  /** */
  Vector<double> currentEvalPt() const 
    {
      return this->ptr()->currentEvalPt();
    }

private:
};
}


#endif
