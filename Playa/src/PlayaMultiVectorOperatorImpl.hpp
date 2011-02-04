/* @HEADER@ */

 /* @HEADER@ */

#ifndef PLAYA_MULTI_VECTOR_OPERATOR_IMPL_HPP
#define PLAYA_MULTI_VECTOR_OPERATOR_IMPL_HPP

#include "PlayaMultiVectorOperatorDecl.hpp"
#include "PlayaVectorImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaSimplifiedLinearOpBaseImpl.hpp"
#endif

namespace Playa
{

/*
 * Construct from an array of vectors and a specifier for the 
 * domain space. 
 */
template <class Scalar> inline
MultiVectorOperator<Scalar>
::MultiVectorOperator(const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain)
  : LinearOpWithSpaces<Scalar>(domain, cols[0].space()),
    cols_(cols)
{
  TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
    "empty multivector given to MultiVectorOperator ctor");
  for (int i=1; i<cols.size(); i++)
  {
    TEST_FOR_EXCEPTION(cols[i].space() != cols[0].space(), std::runtime_error,
      "inconsistent vector spaces in  MultiVectorOperator ctor");
  }
}


/*
 * Apply does an element-by-element multiply between the input 
 * vector, x, and the diagonal values.
 */
template <class Scalar> inline
void MultiVectorOperator<Scalar>
::apply(
  Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in, 
  Vector<Scalar> out) const
{
  if (transApplyType == NO_TRANS)
  {
    out.zero();

    for (int i=0; i<cols_.size(); i++)
    {
      out.update(in[i], cols_[i]);
    }
  }
  else
  {
    out.zero();

    for (int i=0; i<cols_.size(); i++)
    {
      out[i] = in.dot(cols_[i]);
    }
  }
}



/* Return the kth row  */
template <class Scalar> inline
void MultiVectorOperator<Scalar>
::getRow(const int& k, 
  Teuchos::Array<int>& indices, 
  Teuchos::Array<Scalar>& values) const
{
  int low = this->range()->baseGlobalNaturalIndex();
  indices.resize(cols_.size());
  values.resize(cols_.size());
  for (int j=0; j<cols_.size(); j++)
  {
    indices[j] = j;
    values[j] = cols_[j][k-low];
  }
}


  
template <class Scalar> inline
LinearOperator<Scalar> multiVectorOperator(
  const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain)
{
  RCP<LinearOperatorBase<Scalar> > A
    = rcp(new MultiVectorOperator<Scalar>(cols, domain));

  return A;
}


}

#endif
