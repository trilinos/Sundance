/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_ADDED_OP_IMPL_HPP
#define PLAYA_SIMPLE_ADDED_OP_IMPL_HPP



#include "PlayaSimpleAddedOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_Array.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleZeroOpImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;



/*
 * Represent a sum of operators A_0 + A_1 + ... + A_n.
 */
template <class Scalar> inline
SimpleAddedOp<Scalar>::SimpleAddedOp(
  const Array<LinearOperator<Scalar> >& ops)
  : LinearOpWithSpaces<Scalar>(
    ops[0].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (int i=1; i<ops_.size(); i++)
  {
    TEST_FOR_EXCEPT(!(ops[i].range() == ops[0].range()));
    TEST_FOR_EXCEPT(!(ops[i].domain() == ops[0].domain()));
  }
}
  
/* */
template <class Scalar> inline
void SimpleAddedOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleAddedOp::apply()");

  Vector<Scalar> tmp=out.copy();
  tmp.zero();
  for (int i=0; i<ops_.size(); i++)
  {
    Tabs tab1;
    Out::os() << tab1 << "applying term i=" << i << " of " 
              << ops_.size() << std::endl;
    if (transApplyType == Teuchos::NO_TRANS)
      tmp = tmp + ops_[i] * in;
    else if (transApplyType == Teuchos::TRANS)
      tmp = tmp + ops_[i].transpose() * in;
    else 
      TEST_FOR_EXCEPT(transApplyType != Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);
  }
  out.acceptCopyOf(tmp);

  PLAYA_MSG2(this->verb(), tab << "done SimpleAddedOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleAddedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (int i=0; i<ops_.size(); i++)
  {
    if (i > 0) rtn += "+";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
}



template <class Scalar> inline
LinearOperator<Scalar> addedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  /* We will strip out any zero operators */
  Array<LinearOperator<Scalar> > strippedOps;

  for (int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    /* Ignore any zero operators */
    const SimpleZeroOp<Scalar>* zPtr 
      = dynamic_cast<const SimpleZeroOp<Scalar>*>(op_i.ptr().get());

    if (zPtr != 0) continue;

    strippedOps.append(op_i);
  }
  
  TEST_FOR_EXCEPT(strippedOps.size() < 1);
  if (strippedOps.size()==1) return strippedOps[0];
  
  RCP<LinearOperatorBase<Scalar> > op 
    = rcp(new SimpleAddedOp<Scalar>(strippedOps));
  
  return op;
}

template <class Scalar> inline
LinearOperator<Scalar> operator+(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B)
{
  return addedOperator(Array<LinearOperator<Scalar> >(tuple(A, B)));
}

}

#endif
