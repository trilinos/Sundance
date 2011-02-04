/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_COMPOSED_OP_IMPL_HPP
#define PLAYA_SIMPLE_COMPOSED_OP_IMPL_HPP



#include "PlayaSimpleComposedOpDecl.hpp"
#include "PlayaSimpleIdentityOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleIdentityOpImpl.hpp"
#include "PlayaSimpleZeroOpImpl.hpp"
#include "PlayaLinearOpWithSpacesImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;




/*
 * ------ composed operator  
 */
template <class Scalar> inline
SimpleComposedOp<Scalar>::SimpleComposedOp(const Array<LinearOperator<Scalar> >& ops)
  : LinearOpWithSpaces<Scalar>(
    ops[ops.size()-1].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (int i=1; i<ops_.size(); i++)
  {
    TEST_FOR_EXCEPT(!(ops[i].range() == ops[i-1].domain()));
  }
}
  


template <class Scalar> inline
void SimpleComposedOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleComposedOp::apply()");
  if (transApplyType == Teuchos::NO_TRANS)
  {
    Vector<Scalar> tmpIn = in.copy();
    for (int i=0; i<ops_.size(); i++)
    {
      Tabs tab1;
      Vector<Scalar> tmpOut;
      int j = ops_.size()-1-i;
      PLAYA_MSG2(this->verb(), tab1 << "applying op #" << j 
        << " of " << ops_.size());
      ops_[j].apply(tmpIn, tmpOut);
      tmpIn = tmpOut;
    }
    out.acceptCopyOf(tmpIn);
  }

  else if (transApplyType == Teuchos::TRANS)
  {
    Vector<Scalar> tmpIn = in.copy();
    for (int i=0; i<ops_.size(); i++)
    {
      Tabs tab1;
      Vector<Scalar> tmpOut;
      PLAYA_MSG2(this->verb(), tab1 << "applying transpose op #" << i
        << " of " << ops_.size());
      ops_[i].applyTranspose(tmpIn, tmpOut);
      tmpIn = tmpOut;
    }
    out.acceptCopyOf(tmpIn);
  }
  else
  {
    TEST_FOR_EXCEPT(transApplyType != Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);
  }
  PLAYA_MSG2(this->verb(), tab << "done SimpleComposedOp::apply()");
}
  



template <class Scalar> inline
std::string SimpleComposedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (int i=0; i<ops_.size(); i++)
  {
    if (i > 0) rtn += "*";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
}


template <class Scalar> inline
void SimpleComposedOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "ComposedOperator[" << std::endl;
  for (int i=0; i<ops_.size(); i++)
  {
    Tabs tab1;
    os << tab1 << "factor #" << i << std::endl;
    Tabs tab2;
    os << tab2 << ops_[i].description() << std::endl;
    os << std::endl;
  }
  os << tab << "]" <<  std::endl;
}


template <class Scalar> inline
LinearOperator<Scalar> composedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  /* We will strip out any identity operators, and if we find a zero
  * operator the whole works becomes a zero operator */ 
  Array<LinearOperator<Scalar> > strippedOps;

  for (int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    /* if a factor is zero, the whole operator is
     * a zero operator */
    const SimpleZeroOp<Scalar>* zPtr 
      = dynamic_cast<const SimpleZeroOp<Scalar>*>(op_i.ptr().get());

    if (zPtr != 0) 
    {
      VectorSpace<Scalar> r = ops[0].range();
      VectorSpace<Scalar> d = ops[ops.size()-1].domain();
      return zeroOperator(d, r);
    }

    /* if a factor is the identity, skip it */
    const SimpleIdentityOp<Scalar>* IPtr 
      = dynamic_cast<const SimpleIdentityOp<Scalar>*>(op_i.ptr().get());  
    if (IPtr != 0) 
    {
      continue;
    }

    strippedOps.append(op_i);
  }
  
  TEST_FOR_EXCEPT(strippedOps.size() < 1);
  if (strippedOps.size()==1) return strippedOps[0];
  
  RCP<LinearOperatorBase<Scalar> > op 
    = rcp(new SimpleComposedOp<Scalar>(strippedOps));
  return op;
}

template <class Scalar> inline
LinearOperator<Scalar> operator*(const LinearOperator<Scalar>& A, 
  const LinearOperator<Scalar>& B)
{
  return composedOperator(Array<LinearOperator<Scalar> >(tuple(A,B)));
}
  

}

#endif
