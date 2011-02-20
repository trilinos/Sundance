/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_SCALED_OP_IMPL_HPP
#define PLAYA_SIMPLE_SCALED_OP_IMPL_HPP



#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;





/*
 * --- scaled op
 */

template <class Scalar> inline
SimpleScaledOp<Scalar>::SimpleScaledOp(const Scalar& alpha,
  const LinearOperator<Scalar>& A)
  : LinearOpWithSpaces<Scalar>(
    A.domain(), A.range()
    ) 
  , alpha_(alpha), A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleScaledOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleScaledOp::apply()");

  if (transApplyType == Teuchos::NO_TRANS)
    A_.apply(in, out);
  else if (transApplyType == Teuchos::TRANS)
    A_.applyTranspose(in, out);
  else 
    TEST_FOR_EXCEPT(transApplyType !=Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);

  out.scale(alpha_);

  PLAYA_MSG2(this->verb(), tab << "done SimpleScaledOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleScaledOp<Scalar>::description() const 
{
  return "ScaledOp[alpha="  + Teuchos::toString(alpha_)
    + ", " + A_.description() + "]";
}


/* */
template <class Scalar> inline
void SimpleScaledOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "ScaledOp[" << std::endl;
  Tabs tab1;
  os << tab1 << "scale = " << alpha_ << std::endl;
  os << tab1 << "operator = " << A_.description() << std::endl;
  os << tab << "]" << std::endl;
}



template <class Scalar> inline
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RCP<LinearOperatorBase<Scalar> > A 
    = rcp(new SimpleScaledOp<Scalar>(scale, op));

  return A;
}


template <class Scalar> inline
LinearOperator<Scalar> operator*(const Scalar& a, const LinearOperator<Scalar>& A)
{
  return scaledOperator(a, A);
}
  

}

#endif
