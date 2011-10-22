/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_DEFAULT_BLOCK_VECTOR_IMPL_HPP
#define PLAYA_DEFAULT_BLOCK_VECTOR_IMPL_HPP

#include "PlayaDefaultBlockVectorDecl.hpp"
#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaExceptions.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaBlockIteratorImpl.hpp"
#include "PlayaVectorSpaceBaseImpl.hpp"
#include "PlayaDefaultBlockVectorSpaceImpl.hpp"
#include "PlayaBlockVectorBaseImpl.hpp"
#endif





namespace Playa
{
template <class Scalar> inline
DefaultBlockVector<Scalar>
::DefaultBlockVector(const VectorSpace<Scalar>& space)
  : BlockVectorBase<Scalar>(), space_(space), blocks_(space.numBlocks())
{
  for (int i=0; i<space.numBlocks(); i++)
  {
    blocks_[i] = space.getBlock(i).createMember();
  }
}

template <class Scalar> inline
DefaultBlockVector<Scalar>
::DefaultBlockVector(const VectorSpace<Scalar>& space, 
    const Array<Vector<Scalar> >& blocks)
  : BlockVectorBase<Scalar>(), space_(space), blocks_(blocks)
{}

template <class Scalar> inline
void DefaultBlockVector<Scalar>::setBlock(int b, const Vector<Scalar>& block)
{
  PLAYA_BOUNDSCHECK(b, 0, space_.numBlocks(), 
    "DefaultBlockVector::setBlock()");

  TEUCHOS_TEST_FOR_EXCEPTION(
    block.space() != space_.getBlock(b),
    RuntimeError,
    "inconsistent block spaces in setBlock: \n" 
    "old=" << space_.getBlock(b) << std::endl
    << "new=" << block.space());
  
  blocks_[b] = block;
}


template <class Scalar> inline
const Vector<Scalar>& DefaultBlockVector<Scalar>::getBlock(int b) const 
{
  PLAYA_BOUNDSCHECK(b, 0, space_.numBlocks(), 
    "const DefaultBlockVector::setBlock()");

  return blocks_[b];
}

template <class Scalar> inline
Vector<Scalar> DefaultBlockVector<Scalar>::getNonConstBlock(int b) 
{
  PLAYA_BOUNDSCHECK(b, 0, space_.numBlocks(), 
    "non-const DefaultBlockVector::setBlock()");

  return blocks_[b];
}



/** \relates Vector */
template <class Scalar> inline
Vector<Scalar> blockVector(
  const Vector<Scalar>& v1)
{
  Array<Vector<Scalar> > x(1);
  x[0] = v1;
  return blockVector<Scalar>(x);
}


/** \relates Vector */
template <class Scalar> inline
Vector<Scalar> blockVector(
  const Vector<Scalar>& v1,
  const Vector<Scalar>& v2)
{
  Array<Vector<Scalar> > x(2);
  x[0] = v1;
  x[1] = v2;
  return blockVector<Scalar>(x);
}

/** \relates Vector */
template <class Scalar> inline
Vector<Scalar> blockVector(
  const Vector<Scalar>& v1,
  const Vector<Scalar>& v2,
  const Vector<Scalar>& v3)
{
  Array<Vector<Scalar> > x(3);
  x[0] = v1;
  x[1] = v2;
  x[2] = v3;
  return blockVector<Scalar>(x);
}

/** \relates Vector */
template <class Scalar> inline
Vector<Scalar> blockVector(
  const Vector<Scalar>& v1,
  const Vector<Scalar>& v2,
  const Vector<Scalar>& v3,
  const Vector<Scalar>& v4)
{
  Array<Vector<Scalar> > x(4);
  x[0] = v1;
  x[1] = v2;
  x[2] = v3;
  x[3] = v4;
  return blockVector<Scalar>(x);
}

/** \relates Vector */
template <class Scalar> inline
Vector<Scalar> blockVector(const Array<Vector<Scalar> >& x)
{
  Array<VectorSpace<Scalar> > spaces(x.size());
  for (int i=0; i<x.size(); i++) spaces[i] = x[i].space();

  VectorSpace<Scalar> bs = blockSpace(spaces);
  RCP<VectorBase<Scalar> > rtn = rcp(new DefaultBlockVector<Scalar>(bs, x));
  return rtn;
}


  

}

#endif
