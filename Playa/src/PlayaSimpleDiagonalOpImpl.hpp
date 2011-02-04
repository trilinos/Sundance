/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_DIAGONAL_OP_IMPL_HPP
#define PLAYA_SIMPLE_DIAGONAL_OP_IMPL_HPP



#include "PlayaSimpleDiagonalOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOpWithSpacesImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;




/*
 * --- scaled op
 */

template <class Scalar> inline
SimpleDiagonalOp<Scalar>::SimpleDiagonalOp(
  const Vector<Scalar>& diag)
  : LinearOpWithSpaces<Scalar>(
    diag.space(), diag.space()
    ), diag_(diag)
{}
  
/* */
template <class Scalar> inline
void SimpleDiagonalOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleDiagonalOp::apply()");

  Vector<Scalar> tmp = in.dotStar(diag_);
  out.acceptCopyOf(tmp);

  PLAYA_MSG2(this->verb(), tab << "done SimpleDiagonalOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleDiagonalOp<Scalar>::description() const 
{
  return "DiagonalOp[diag=" + diag_.description() + "]";
}


/* */
template <class Scalar> inline
void SimpleDiagonalOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "DiagonalOp[" << std::endl;
  Tabs tab1;
  os << tab1 << "diag = " << diag_ << std::endl;
  os << tab << "]" << std::endl;
}



template <class Scalar> inline
LinearOperator<Scalar> diagonalOperator(
  const Vector<Scalar>& diag)
{
  RCP<LinearOperatorBase<Scalar> > A 
    = rcp(new SimpleDiagonalOp<Scalar>(diag));

  return A;
}



}

#endif
