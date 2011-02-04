/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_TRANSPOSED_OP_IMPL_HPP
#define PLAYA_SIMPLE_TRANSPOSED_OP_IMPL_HPP



#include "PlayaSimpleTransposedOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleZeroOpImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;




/*
 * --- transposed op
 */

template <class Scalar> inline
SimpleTransposedOp<Scalar>::SimpleTransposedOp(const LinearOperator<Scalar>& A)
  : LinearOpWithSpaces<Scalar>(
    A.range(), A.domain()
    ) 
  , A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleTransposedOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleTransposedOp::apply()");

  if (transApplyType == Teuchos::NO_TRANS)
    A_.applyTranspose(in, out);
  else if (transApplyType == Teuchos::TRANS)
    A_.apply(in, out);
  else 
    TEST_FOR_EXCEPT(transApplyType !=Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);

  PLAYA_MSG2(this->verb(), tab << "done SimpleTransposedOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleTransposedOp<Scalar>::description() const 
{
  return "(" + A_.description() + "^T)";
}



template <class Scalar> inline
LinearOperator<Scalar> transposedOperator(
  const LinearOperator<Scalar>& op)
{

  /* If the operator is a transpose, return the untransposed op */
  const SimpleTransposedOp<Scalar>* tPtr
    = dynamic_cast<const SimpleTransposedOp<Scalar>*>(op.ptr().get());
  if (tPtr)
  {
    return tPtr->op();
  }

  /* If the operator is zero, return a transposed zero */
  const SimpleZeroOp<Scalar>* zPtr 
    = dynamic_cast<const SimpleZeroOp<Scalar>*>(op.ptr().get());

  if (zPtr != 0) 
  {
    VectorSpace<Scalar> r = op.range();
    VectorSpace<Scalar> d = op.domain();
    return zeroOperator(r, d);
  }


  /* Return a transposed operator */
  RCP<LinearOperatorBase<Scalar> > A
    = rcp(new SimpleTransposedOp<Scalar>(op));
      
  return A;
}

  


}

#endif
