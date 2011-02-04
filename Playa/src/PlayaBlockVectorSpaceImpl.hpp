/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BLOCKEDVECTORSPACEBASEIMPL_HPP
#define PLAYA_BLOCKEDVECTORSPACEBASEIMPL_HPP

#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"

namespace Playa
{
template <class Scalar> inline
int BlockVectorSpaceBase<Scalar>::dim() const 
{
  int rtn = 0;
  for (int i=0; i<this->numBlocks(); i++) 
  {
    rtn += this->getBlock(i).dim();
  }
  return rtn;
}

template <class Scalar> inline
int BlockVectorSpaceBase<Scalar>::numLocalElements() const 
{
  int rtn = 0;
  for (int i=0; i<this->numBlocks(); i++) 
  {
    rtn += this->getBlock(i).numLocalElements();
  }
  return rtn;
}

template <class Scalar> inline
bool BlockVectorSpaceBase<Scalar>::isCompatible(const VectorSpaceBase<Scalar>* other) const 
{
  if (this->numBlocks() != other->numBlocks()) return false;

  const BlockVectorSpaceBase<Scalar>* bs 
    = dynamic_cast<const BlockVectorSpaceBase<Scalar>* >(other);
    
  if (bs == 0) return false;

  for (int i=0; i<this->numBlocks(); i++) 
  {
    if (! (this->getBlock(i).isCompatible(bs->getBlock(i))))
      return false;
  }
  return true;
}


template <class Scalar> inline
std::string BlockVectorSpaceBase<Scalar>::description() const 
{
  std::ostringstream rtn;
  rtn << "BlockVS[";
  for (int b=0; b<this->numBlocks(); b++)
  {
    if (b > 0) rtn << ", ";
    rtn << this->getBlock(b).description();
  }
  rtn << "]";
  return rtn.str();
}


  
  
}

#endif
