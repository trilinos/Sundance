/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORSPACEIMPL_HPP
#define PLAYA_VECTORSPACEIMPL_HPP

#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaBlockVectorSpaceImpl.hpp"
#endif

using namespace Teuchos;

namespace Playa
{


static inline Time& createVecTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("vector allocation"); 
  return *rtn;
}

 
//========================================================================
template <class Scalar> inline
bool VectorSpace<Scalar>::operator==(const VectorSpace<Scalar>& other) const 
{
  return isCompatible(other);  
}


//========================================================================
template <class Scalar> inline
bool VectorSpace<Scalar>::operator!=(const VectorSpace<Scalar>& other) const 
{
  return !(operator==(other));
}
    


//========================================================================
template <class Scalar> inline
Vector<Scalar> VectorSpace<Scalar>::createMember() const 
{
  Vector<Scalar> rtn = this->ptr()->createMember(*this);
  rtn.setToConstant(0.0);
  return rtn;
}
    


//========================================================================
template <class Scalar> inline
int VectorSpace<Scalar>::baseGlobalNaturalIndex() const
{
  return this->ptr()->baseGlobalNaturalIndex();
}

//========================================================================
template <class Scalar> inline
int VectorSpace<Scalar>::numLocalElements() const
{
  return this->ptr()->numLocalElements();
}
    



//========================================================================
template <class Scalar> inline
bool VectorSpace<Scalar>::isCompatible(const VectorSpace<Scalar>& vecSpc) const 
{
  TEST_FOR_EXCEPTION(vecSpc.ptr().get() == 0, std::runtime_error,
                     "null argument in VectorSpace<Scalar>::isCompatible()");
  return this->ptr().get()->isCompatible(vecSpc.ptr().get());
}





//========================================================================
template <class Scalar> inline
bool VectorSpace<Scalar>::contains(const Vector<Scalar> &vec) const
{
  return (operator==(vec.space()));
}


//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::numBlocks() const
{
  return this->ptr()->numBlocks();
}


//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::isBlockSpace() const
{
  return dynamic_cast<const BlockVectorSpaceBase<Scalar>*>(this->ptr().get())!=0;
}



//========================================================================
template <class Scalar> inline
const VectorSpace<Scalar>& VectorSpace<Scalar>::getBlock(const int i) const
{
  const BlockVectorSpaceBase<Scalar>* bvs = 
    dynamic_cast<const BlockVectorSpaceBase<Scalar>* > (this->ptr().get());
  TEST_FOR_EXCEPTION(bvs == 0 && numBlocks()!=1, std::runtime_error,
    "getBlock called for a space that "
    "is not a BlockVectorSpace" << std::endl);
  if (bvs != 0)
    {
      return bvs->getBlock(i);
    }
  return *this;
}

//========================================================================
template <class Scalar> inline
const VectorSpace<Scalar>& VectorSpace<Scalar>
::getBlock(const BlockIterator<Scalar>& b) const
{
  /* Check that the block iterator is valid */
  TEST_FOR_EXCEPTION(b.atEnd(), RuntimeError, 
    "Attempt to use a block iterator that's run off end");

  return this->getBlock(b.blockIndex());
}


//========================================================================
template <class Scalar> inline
const VectorSpace<Scalar>& VectorSpace<Scalar>
::getBlock(const std::deque<int>& b) const
{
  /* Check that the block iterator is valid */
  if (b.size()==0) return *this;
  
  if (b.size()==1) 
  {
    return this->getBlock(b.front());
  }

  int b0 = b.front();
  std::deque<int> tmp = b;
  tmp.pop_front();
  return this->getBlock(b0).getBlock(tmp);
}





template <class Scalar> inline
BlockIterator<Scalar> VectorSpace<Scalar>::beginBlock() const
{
  return BlockIterator<Scalar>(*this, false);
}

template <class Scalar> inline
BlockIterator<Scalar> VectorSpace<Scalar>::endBlock() const
{
  return BlockIterator<Scalar>(*this, true);
}


template <class Scalar> inline
int VectorSpace<Scalar>::mapToGNI(
  const BlockIterator<Scalar>& b,
  int indexWithinBlock) const 
{
  int rtn = this->baseGlobalNaturalIndex();
  if (this->numBlocks() > 1)
  {
    for (BlockIterator<Scalar> a=this->beginBlock(); a < b; a++)
    {
      rtn += getBlock(a).numLocalElements(); 
    }
  }
  return rtn + indexWithinBlock;
}

template <class Scalar> inline
bool VectorSpace<Scalar>::containsGNI(int gni) const
{
  return (gni >= this->baseGlobalNaturalIndex()
    && gni < this->baseGlobalNaturalIndex()+this->numLocalElements());
}

template <class Scalar> inline
void VectorSpace<Scalar>::getBlockAndOffsetFromGNI(int gni,
  BlockIterator<Scalar>& block, int& indexWithinBlock) const 
{
  int low = this->baseGlobalNaturalIndex();
  
  if (this->numBlocks() == 1)
  {
    indexWithinBlock = gni - low;
  }
  else
  {
    int sizeSum = low;
    for (BlockIterator<Scalar> b=this->beginBlock(); 
         b != this->endBlock(); b++)
    {
      if (this->getBlock(b).containsGNI(gni))
      {
        block = b;
        indexWithinBlock = gni - sizeSum;
      }
      sizeSum += this->getBlock(b).numLocalElements();
    }
  }
}

}



#endif
