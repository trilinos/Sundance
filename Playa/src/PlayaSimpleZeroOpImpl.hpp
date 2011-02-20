/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_ZERO_OP_IMPL_HPP
#define PLAYA_SIMPLE_ZERO_OP_IMPL_HPP



#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;





/* ---- Zero op ------- */

template <class Scalar> inline
SimpleZeroOp<Scalar>::SimpleZeroOp(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : LinearOpWithSpaces<Scalar>(domain, range) {}



template <class Scalar> inline
void SimpleZeroOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleZeroOp::apply()");

  out.zero();

  PLAYA_MSG2(this->verb(), tab << "done SimpleZeroOp::apply()");
}

/* */
template <class Scalar> inline
std::string SimpleZeroOp<Scalar>::description() const 
{return "ZeroOp(domain=" 
    + this->domain()->description() 
    + ", range=" + this->range()->description() + ")";}



template <class Scalar> inline
LinearOperator<Scalar> zeroOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
{
  RCP<LinearOperatorBase<Scalar> > op 
    = rcp(new SimpleZeroOp<Scalar>(domain, range));

  return op;
}


}

#endif
