/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_IDENTITY_OP_IMPL_HPP
#define PLAYA_SIMPLE_IDENTITY_OP_IMPL_HPP



#include "PlayaSimpleIdentityOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearOpWithSpacesImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;



/* ---- Identity op ------- */

template <class Scalar> inline
SimpleIdentityOp<Scalar>::SimpleIdentityOp(const VectorSpace<Scalar>& space)
  : LinearOpWithSpaces<Scalar>(space, space) {}


template <class Scalar> inline
void SimpleIdentityOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleZeroOp::apply()");
  out.acceptCopyOf(in);
  PLAYA_MSG2(this->verb(), tab << "done SimpleIdentityOp::apply()");
}

template <class Scalar> inline
std::string SimpleIdentityOp<Scalar>::description() const 
{return "I(" + this->domain()->description() + ")";}



template <class Scalar> inline
LinearOperator<Scalar> identityOperator(
  const VectorSpace<Scalar>& space)
{
  RCP<LinearOperatorBase<Scalar> > op 
    = rcp(new SimpleIdentityOp<Scalar>(space));

  return op;
}



}

#endif
