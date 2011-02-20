/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_DEFAULT_BLOCK_VECTOR_IMPL_HPP
#define PLAYA_DEFAULT_BLOCK_VECTOR_IMPL_HPP

#include "PlayaDefaultBlockVectorDecl.hpp"
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
void DefaultBlockVector<Scalar>::setBlock(int b, const Vector<Scalar>& block)
{
  PLAYA_BOUNDSCHECK(b, 0, space_.numBlocks(), 
    "DefaultBlockVector::setBlock()");

  TEST_FOR_EXCEPTION(
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



  

}

#endif
