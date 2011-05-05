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
void BlockVectorBase<Scalar>::update(const Scalar& alpha, 
  const VectorBase<Scalar>* other,
  const Scalar& gamma)
{
  const BlockVectorBase<Scalar>* bvo 
    = dynamic_cast<const BlockVectorBase<Scalar>* >(other);
  for (int b=0; b<this->numBlocks(); b++)
  {
    this->getNonConstBlock(b).update(alpha, bvo->getBlock(b), gamma);
  }
}

template <class Scalar> inline
void BlockVectorBase<Scalar>::update(
  const Scalar& alpha, const VectorBase<Scalar>* x,
  const Scalar& beta, const VectorBase<Scalar>* y,
  const Scalar& gamma)
{
  const BlockVectorBase<Scalar>* bx
    = dynamic_cast<const BlockVectorBase<Scalar>* >(x);
  const BlockVectorBase<Scalar>* by 
    = dynamic_cast<const BlockVectorBase<Scalar>* >(y);

  for (int b=0; b<this->numBlocks(); b++)
  {
    this->getNonConstBlock(b).ptr()->update(
      alpha, bx->getBlock(b).ptr().get(), 
      beta, by->getBlock(b).ptr().get(),
      gamma);
  }
}


template <class Scalar> inline
void BlockVectorBase<Scalar>::update(
  const Scalar& alpha, const VectorBase<Scalar>* x,
  const Scalar& beta, const VectorBase<Scalar>* y,
  const Scalar& gamma, const VectorBase<Scalar>* z,
  const Scalar& delta)
{
  const BlockVectorBase<Scalar>* bx
    = dynamic_cast<const BlockVectorBase<Scalar>* >(x);
  const BlockVectorBase<Scalar>* by 
    = dynamic_cast<const BlockVectorBase<Scalar>* >(y);
  const BlockVectorBase<Scalar>* bz 
    = dynamic_cast<const BlockVectorBase<Scalar>* >(z);

  for (int b=0; b<this->numBlocks(); b++)
  {
    this->getNonConstBlock(b).ptr()->update(
      alpha, bx->getBlock(b).ptr().get(), 
      beta, by->getBlock(b).ptr().get(),
      gamma, bz->getBlock(b).ptr().get(),
      delta);
  }
}

template <class Scalar> inline
Scalar BlockVectorBase<Scalar>::dot(
  const VectorBase<Scalar>* other) const 
{
  const BlockVectorBase<Scalar>* bx
    = dynamic_cast<const BlockVectorBase<Scalar>* >(other);

  Scalar rtn = 0.0;

  for (int b=0; b<this->numBlocks(); b++)
  {
    rtn += this->getBlock(b).ptr()->dot(bx->getBlock(b).ptr().get());
  }
  return rtn;
}

template <class Scalar> inline
Scalar BlockVectorBase<Scalar>::norm2() const 
{
  Scalar rtn = 0.0;

  for (int b=0; b<this->numBlocks(); b++)
  {
    Scalar tmp = 0.0;
    tmp = this->getBlock(b).ptr()->norm2();
    rtn += tmp*tmp;
  }
  return sqrt(rtn);
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
