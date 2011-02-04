/* @HEADER@ */
//
/* @HEADER@ */

#ifndef PLAYA_INVERSEOPERATOR_IMPL_HPP
#define PLAYA_INVERSEOPERATOR_IMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaTabs.hpp"
#include "PlayaSolverState.hpp"
#include "PlayaInverseOperatorDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"

#include "Teuchos_RefCountPtr.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#endif


namespace Playa
{
using Teuchos::RCP;

/*
 * Ctor with a linear operator and a solver specified.
 */
template <class Scalar> inline
InverseOperator<Scalar>::InverseOperator(const LinearOperator<Scalar>& op, 
  const LinearSolver<Scalar>& solver)
  : LinearOpWithSpaces<Scalar>(op.domain(), op.range()), 
    op_(op), solver_(solver) {;}


template <class Scalar> inline
void InverseOperator<Scalar>::apply(
  Teuchos::ETransp applyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const 
{
  
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "InverseOperator::generalApply()");

  TEST_FOR_EXCEPTION(dynamic_cast<SimpleZeroOp<Scalar>* >(op_.ptr().get()) != 0, std::runtime_error,
    "InverseOperator<Scalar>::apply() called on a zero operator.");

  TEST_FOR_EXCEPTION(op_.domain().dim() != op_.range().dim(), std::runtime_error,
    "InverseOperator<Scalar>::apply() called on a non-square operator.");

  Vector<Scalar> result = out.space().createMember();

  SolverState<Scalar> haveSoln;
  if (applyType==NO_TRANS)
  {
    haveSoln = solver_.solve(op_, in, result);
  }
  else
  {
    haveSoln = solver_.solve(op_.transpose(), in, result);
  }


  TEST_FOR_EXCEPTION(haveSoln.finalState() != SolveConverged, 
    std::runtime_error,
    "InverseOperator<Scalar>::apply() " 
    << haveSoln.stateDescription());

  out.acceptCopyOf(result);

  PLAYA_MSG2(this->verb(), tab << "done InverseOperator::generalApply()");
}



template <class Scalar> 
void InverseOperator<Scalar>::print(std::ostream& os) const
{
  Tabs tab(0);
  os << tab << "InverseOperator[" << std::endl;
  Tabs tab1;
  os << tab1 << "op=" << op_ << std::endl;
  os << tab << "]" << std::endl;
}


template <class Scalar> 
LinearOperator<Scalar> 
inverse(const LinearOperator<Scalar>& op, 
  const LinearSolver<Scalar>& solver)
{
  /* check for the case (A^-1)^-1 */
  const InverseOperator<Scalar>* invOp 
    = dynamic_cast<const InverseOperator<Scalar>*>(op.ptr().get());
  if (invOp) return invOp->op();

  /* otherwise return an implicit inverse operator */
  RCP<LinearOperatorBase<Scalar> > rtn 
    = rcp(new InverseOperator<Scalar>(op, solver));
  return rtn;
}

}

#endif
