/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BLOCK_VECTOR_BASE_IMPL_HPP
#define PLAYA_BLOCK_VECTOR_BASE_IMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaTabs.hpp"
#include "PlayaBlockVectorBaseDecl.hpp"
#include "PlayaVectorDecl.hpp"

namespace Playa
{

template <class Scalar> inline 
bool BlockVectorBase<Scalar>::hasMoreChunks() const 
{
  if (currentBlock_ < this->numBlocks()-1)
  {
    return true;
  } 
  if (currentBlock_ == this->numBlocks()-1)
  {
    return this->getBlock(currentBlock_).hasMoreChunks();
  }
  return false;
}

template <class Scalar> inline 
ConstDataChunk<Scalar> BlockVectorBase<Scalar>::nextConstChunk() const
{
  if (!this->getBlock(currentBlock_).hasMoreChunks())
  {
    currentBlock_++;
  }
  return this->getBlock(currentBlock_).nextConstChunk();
}


template <class Scalar> inline 
NonConstDataChunk<Scalar> BlockVectorBase<Scalar>::nextChunk()
{
  if (!this->getBlock(currentBlock_).hasMoreChunks())
  {
    currentBlock_++;
  }
  return this->getNonConstBlock(currentBlock_).nextChunk();
}

template <class Scalar> inline
void BlockVectorBase<Scalar>::rewind() const
{
  Tabs tab;
  currentBlock_ = 0;
  for (int b=0; b<this->numBlocks(); b++)
  {
    this->getBlock(b).rewind();
  }
}

template <class Scalar> inline
std::string BlockVectorBase<Scalar>::description() const
{
  std::ostringstream oss;
  oss << "BlockVector[";
  for (int i=0; i<this->numBlocks(); i++) 
  {
    if (i > 0) oss << ", ";
    oss << this->getBlock(i).description();
  }
  oss << "]";
  return oss.str();
}

  
}

#endif
