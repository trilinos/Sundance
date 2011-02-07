/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORIMPL_HPP
#define PLAYA_VECTORIMPL_HPP


#include "PlayaVectorDecl.hpp"
#include "PlayaBlockVectorBaseDecl.hpp"
#include "PlayaVectorFunctorsImpl.hpp"
#include "PlayaSingleChunkVector.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaPrintable.hpp"
#include "PlayaExceptions.hpp"
#include "PlayaGeneralizedIndex.hpp"
#include "PlayaOut.hpp"
#include "Teuchos_StrUtils.hpp"
#include "PlayaTabs.hpp"
#include "PlayaDebug.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorSpaceImpl.hpp"
#endif




namespace Playa
{
using Playa::Out;
using Playa::Tabs;
using Playa::Printable;

//===========================================================================
template <class Scalar> 
int Vector<Scalar>::numBlocks() const 
{
  return this->ptr()->numBlocks();
}  


//===========================================================================
template <class Scalar> 
void Vector<Scalar>::setBlock(int i, const Vector<Scalar>& v)
{
  BlockVectorBase<Scalar>* bv = 
    dynamic_cast<BlockVectorBase<Scalar>* >(this->ptr().get());
  TEST_FOR_EXCEPTION(bv == 0, std::runtime_error,
    "setBlock() called on a vector is not a block vector");
  bv->setBlock(i, v);
}  


//===========================================================================
template <class Scalar> 
const Vector<Scalar>& Vector<Scalar>::getBlock(int i) const
{
  const BlockVectorBase<Scalar>* bv = 
    dynamic_cast<const BlockVectorBase<Scalar>* >(this->ptr().get());
  if (bv==0 && numBlocks()==1) return *this;

  TEST_FOR_EXCEPTION(bv == 0, std::runtime_error,
    "getBlock() called on a vector is not a block vector");
  return bv->getBlock(i);
}

//===========================================================================
template <class Scalar> 
Vector<Scalar> Vector<Scalar>::getNonConstBlock(int i)
{
  BlockVectorBase<Scalar>* bv = 
    dynamic_cast<BlockVectorBase<Scalar>* >(this->ptr().get());
  if (bv==0 && numBlocks()==1) return *this;

  TEST_FOR_EXCEPTION(bv == 0, std::runtime_error,
    "getBlock() called on a vector is not a block vector");
  return bv->getNonConstBlock(i);
}


//========================================================================
template <class Scalar> inline
const Vector<Scalar>& Vector<Scalar>
::getBlock(const BlockIterator<Scalar>& b) const
{
  /* Check that the block iterator is valid */
  TEST_FOR_EXCEPTION(b.atEnd(), RuntimeError, 
    "Attempt to use a block iterator that's run off end");

  return this->getBlock(b.blockIndex());
}



//========================================================================
template <class Scalar> inline
Vector<Scalar> Vector<Scalar>
::getNonConstBlock(const BlockIterator<Scalar>& b)
{
  /* Check that the block iterator is valid */
  TEST_FOR_EXCEPTION(b.atEnd(), RuntimeError, 
    "Attempt to use a block iterator that's run off end");

  return this->getNonConstBlock(b.blockIndex());
}


//========================================================================
template <class Scalar> inline
const Vector<Scalar>& Vector<Scalar>
::getBlock(const std::deque<int>& b) const
{
  /* Check that the block iterator is valid */
  TEST_FOR_EXCEPTION(b.size()==0, RuntimeError, 
    "Attempt to use an empty block iterator");
  
  if (b.size()==1) 
  {
    return this->getBlock(b.front());
  }

  int b0 = b.front();
  std::deque<int> tmp = b;
  tmp.pop_front();
  return this->getBlock(b0).getBlock(tmp);
}

//========================================================================
template <class Scalar> inline
Vector<Scalar> Vector<Scalar>
::getNonConstBlock(const std::deque<int>& b) 
{
  /* Check that the block iterator is valid */
  TEST_FOR_EXCEPTION(b.size()==0, RuntimeError, 
    "Attempt to use an empty block iterator");
  
  if (b.size()==1) 
  {
    return this->getNonConstBlock(b.front());
  }

  int b0 = b.front();
  std::deque<int> tmp = b;
  tmp.pop_front();
  return this->getNonConstBlock(b0).getNonConstBlock(tmp);
}

//===========================================================================
template <class Scalar> 
ConstDataChunk<Scalar> Vector<Scalar>::nextConstChunk() const
{
  return this->ptr()->nextConstChunk();
}

//===========================================================================
template <class Scalar> 
NonConstDataChunk<Scalar> Vector<Scalar>::nextChunk() 
{
  return this->ptr()->nextChunk();
}


//===========================================================================
template <class Scalar> 
bool Vector<Scalar>::hasMoreChunks() const 
{
  return this->ptr()->hasMoreChunks();
}

//===========================================================================
template <class Scalar> 
void Vector<Scalar>::rewind() const
{
  return this->ptr()->rewind();
}



//===========================================================================

template <class Scalar> 
void Vector<Scalar>::print(std::ostream& os) const 
{
  const Printable* p = 
    dynamic_cast<const Printable* >(this->ptr().get());
  if (false && p != 0)
  {
    p->print(os);
    return;
  }
  else
  {
    Tabs tab(0);
    int np = this->comm().getNProc();
    int rank = this->comm().getRank();
    if (rank==0) 
    {
      os << tab << this->description() << std::endl;
    }
    this->comm().synchronize();        
    
    os << tab << "rank= " << std::setw(10) << rank << " base GNI=" 
       << std::setw(10) << this->space().baseGlobalNaturalIndex() 
       << " local size=" << std::setw(10) 
       << this->space().numLocalElements() << std::endl; 
    for (BlockIterator<Scalar> b=this->space().beginBlock(); 
         b!=this->space().endBlock(); b++)
    {
      Tabs tab1;
      this->comm().synchronize();        
      if (rank==0 && this->numBlocks()>1)
        os << tab1 << "Block=" << b << std::endl;
      for (int r=0; r<np; r++)
      {
        Tabs tab2;
        this->comm().synchronize();        
        if (rank != r) continue;
        if (np > 1) os << tab2 << "Processor=" << r << std::endl;        
    
        const Vector<Scalar>& xBlock = this->getBlock(b);        
        os << tab2 << "rank= " << rank << " base GNI=" 
           << xBlock.space().baseGlobalNaturalIndex() << std::endl; 

        int low = xBlock.space().baseGlobalNaturalIndex();
        while(xBlock.hasMoreChunks())
        {
          Tabs tab3;
          ConstDataChunk<Scalar> myChunk = xBlock.nextConstChunk();
          const Scalar* me = myChunk.values();
          for (int i=0; i<myChunk.size(); i++)
          {
            os << tab3 << std::setw(10) << low+i << " " << std::setw(20)
               << me[i] << std::endl;
          }
        }
        for(int i=0; i<6; i++)
        {
          this->comm().synchronize();        
        }
      }
    }
    
    this->rewind();
  }
}

template <class Scalar> 
std::string Vector<Scalar>::description() const 
{
  const Describable* d = 
    dynamic_cast<const Describable* >(this->ptr().get());
  if (d != 0)
  {
    return d->description();
  }
  return "Vector[type=unknown, dim=" + Teuchos::toString(dim()) + "]";
}

  

//===========================================================================

template <class Scalar>
template <class UF> inline
Vector<Scalar>& Vector<Scalar>::applyUnaryFunctor(const UF& func) 
{
  TimeMonitor t(*opTimer());

  for (BlockIterator<Scalar> b=this->space().beginBlock(); b!=this->space().endBlock(); b++)
  {
    Vector<Scalar> xBlock = this->getNonConstBlock(b);
    while(xBlock.hasMoreChunks())
    {
      NonConstDataChunk<Scalar> myChunk = xBlock.nextChunk();
      Scalar* me = myChunk.values();
      for (int i=0; i<myChunk.size(); i++)
      {
        me[i] = func(me[i]);
      }
    }
  }

  this->rewind();

  return *this;
}

//===========================================================================

template <class Scalar>
template <class UF> inline
Vector<Scalar>& Vector<Scalar>::acceptUnaryFunctor(const UF& func,
  const Vector<Scalar>& other) 
{
  TimeMonitor t(*opTimer());

  for (BlockIterator<Scalar> b=this->space().beginBlock(); b!=this->space().endBlock(); b++)
  {
    Vector<Scalar> myBlock = this->getNonConstBlock(b);
    Vector<Scalar> yourBlock = other.getBlock(b);
    while(myBlock.hasMoreChunks())
    {
      NonConstDataChunk<Scalar> myChunk = myBlock.nextChunk();
      ConstDataChunk<Scalar> yourChunk = yourBlock.nextConstChunk();
      Scalar* me = myChunk.values();
      const Scalar* you = yourChunk.values();
      for (int i=0; i<myChunk.size(); i++)
      {
        me[i] = func(you[i]);
      }
    }
  }

  other.rewind();
  this->rewind();

  return *this;
}

//===========================================================================

template <class Scalar>
template <class VF> inline
Vector<Scalar>& Vector<Scalar>::applyBinaryFunctor(const VF& func,
  const Vector<Scalar>& other) 
{
  TimeMonitor t(*opTimer());

  for (BlockIterator<Scalar> b=this->space().beginBlock(); b!=this->space().endBlock(); b++)
  {
    Vector<Scalar> myBlock = this->getNonConstBlock(b);
    Vector<Scalar> yourBlock = other.getBlock(b);
    while(myBlock.hasMoreChunks())
    {
      NonConstDataChunk<Scalar> myChunk = myBlock.nextChunk();
      ConstDataChunk<Scalar> yourChunk = yourBlock.nextConstChunk();
      Scalar* me = myChunk.values();
      const Scalar* you = yourChunk.values();
      for (int i=0; i<myChunk.size(); i++)
      {
        me[i] = func(me[i], you[i]);
      }
    }
  }

  other.rewind();
  this->rewind();

  return *this;
}

//===========================================================================
template <class Scalar> 
template <class VF> inline
Vector<Scalar>& Vector<Scalar>::applyTernaryFunctor(
  const VF& func,
  const Vector<Scalar>& y, 
  const Vector<Scalar>& z)
{
  TimeMonitor t(*opTimer());

  for (BlockIterator<Scalar> b=this->space().beginBlock(); b!=this->space().endBlock(); b++)
  {
    Vector<Scalar> xBlock = this->getNonConstBlock(b);
    Vector<Scalar> yBlock = y.getBlock(b);
    Vector<Scalar> zBlock = z.getBlock(b);
    while(xBlock.hasMoreChunks())
    {
      NonConstDataChunk<Scalar> xChunk = xBlock.nextChunk();
      ConstDataChunk<Scalar> yChunk = yBlock.nextConstChunk();
      ConstDataChunk<Scalar> zChunk = zBlock.nextConstChunk();
      Scalar* xv = xChunk.values();
      const Scalar* yv = yChunk.values();
      const Scalar* zv = zChunk.values();
      for (int i=0; i<xChunk.size(); i++)
      {
        xv[i] = func(xv[i], yv[i], zv[i]);
      }
    }
  }

  this->rewind();
  y.rewind();
  z.rewind();

  return *this;
}



//===========================================================================

template <class Scalar>
template <class RF> inline
typename PlayaFunctors::VectorFunctorTraits<Scalar, RF>::ReturnType
Vector<Scalar>::applyUnaryReductionFunctor(const RF& func) const 
{
  TimeMonitor t(*opTimer());

  for (BlockIterator<Scalar> b=this->space().beginBlock(); b!=this->space().endBlock(); b++)
  {
    Vector<Scalar> xBlock = this->getBlock(b);
    while(xBlock.hasMoreChunks())
    {
      ConstDataChunk<Scalar> myChunk = xBlock.nextConstChunk();
      const Scalar* me = myChunk.values();
      for (int i=0; i<myChunk.size(); i++)
      {
        func.step(i, me[i]);
      }
    }
  }

  func.postProc();
  this->rewind();

  return func.result();
}


//===========================================================================

template <class Scalar>
template <class RF> inline
typename PlayaFunctors::VectorFunctorTraits<Scalar, RF>::ReturnType
Vector<Scalar>::applyBinaryReductionFunctor(const RF& func, const Vector<Scalar>& y) const 
{
  TimeMonitor t(*opTimer());

  for (BlockIterator<Scalar> b=this->space().beginBlock(); b!=this->space().endBlock(); b++)
  {
    Vector<Scalar> xBlock = this->getBlock(b);
    Vector<Scalar> yBlock = y.getBlock(b);
    while(xBlock.hasMoreChunks())
    {
      ConstDataChunk<Scalar> xChunk = xBlock.nextConstChunk();
      ConstDataChunk<Scalar> yChunk = yBlock.nextConstChunk();
      const Scalar* me = xChunk.values();
      const Scalar* you = yChunk.values();
      for (int i=0; i<xChunk.size(); i++)
      {
        func.step(i, me[i], you[i]);
      }
    }
  }

  func.postProc();
  this->rewind();
  y.rewind();

  return func.result();
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::scale(const Scalar& alpha)
{
  return applyUnaryFunctor(PlayaFunctors::ScalarMult<Scalar>(alpha));
}

//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::reciprocal()
{
  return applyUnaryFunctor(PlayaFunctors::Reciprocal<Scalar>());
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::abs()
{
  return applyUnaryFunctor(PlayaFunctors::Abs<Scalar>());
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::randomize()
{
  return applyUnaryFunctor(PlayaFunctors::Random<Scalar>());
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
  const Vector<Scalar>& x, const Scalar& gamma)
{
  return applyBinaryFunctor(PlayaFunctors::LC2<Scalar>(gamma, alpha), x);
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
{
  if (this->ptr().get()==0 || !this->space().isCompatible(x.space()) )
  {
    this->ptr() = x.space().createMember().ptr();
  }
  return acceptUnaryFunctor(PlayaFunctors::Identity<Scalar>(), x);
}

template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::copy() const 
{
  Vector<Scalar> rtn = space().createMember();
  rtn.acceptCopyOf(*this);
  return rtn;
}



//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::dotStar(const Vector<Scalar>& other) 
{
  return applyBinaryFunctor(PlayaFunctors::DotStar<Scalar>(), other);
}

//===========================================================================
template <class Scalar> inline 
Vector<Scalar>& Vector<Scalar>::dotSlash(const Vector<Scalar>& other) 
{
  return applyBinaryFunctor(PlayaFunctors::DotSlash<Scalar>(), other);
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::dotStar(const Vector<Scalar>& other) const
{
  Vector<Scalar> rtn = space().createMember();
  rtn.acceptCopyOf(*this);
  rtn.dotStar(other);
  return rtn;
}

//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::dotSlash(const Vector<Scalar>& other) const
{
  Vector<Scalar> rtn = space().createMember();
  rtn.acceptCopyOf(*this);
  rtn.dotSlash(other);
  return rtn;
}


//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::abs() const 
{
  Vector<Scalar> rtn = space().createMember();
  rtn.acceptCopyOf(*this);
  rtn.abs();
  return rtn;
}





//===========================================================================
template <class Scalar> inline 
Vector<Scalar> Vector<Scalar>::reciprocal() const 
{
  Vector<Scalar> rtn = space().createMember();
  rtn.acceptCopyOf(*this);
  rtn.reciprocal();
  return rtn;
}



//===========================================================================
template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, 
  const Vector<Scalar>& x, 
  const Scalar& beta, 
  const Vector<Scalar>& y, 
  const Scalar& gamma)
{
  return applyTernaryFunctor(PlayaFunctors::LC3<Scalar>(gamma, alpha, beta), x, y);
}



//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::dot(const Vector<Scalar>& other) const 
{
  return applyBinaryReductionFunctor(
    PlayaFunctors::DotProduct<Scalar>(this->comm()), other);
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::operator*(const Vector<Scalar>& other) const 
{
  return dot(other);
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm1() const 
{
  return applyUnaryReductionFunctor(PlayaFunctors::Norm1<Scalar>(this->comm()));
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2() const 
{
  return applyUnaryReductionFunctor(PlayaFunctors::Norm2<Scalar>(this->comm()));
}




//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::norm2(const Vector<Scalar>& weights) const 
{
  return applyBinaryReductionFunctor(PlayaFunctors::WeightedNorm2<Scalar>(this->comm()), weights);
}





//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::normInf() const 
{
  return applyUnaryReductionFunctor(PlayaFunctors::NormInf<Scalar>(this->comm()));
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::zero()
{
  setToConstant(0.0);
}




//===========================================================================
template <class Scalar> inline 
void Vector<Scalar>::setToConstant(const Scalar& alpha)
{
  applyUnaryFunctor(PlayaFunctors::SetConstant<Scalar>(alpha));
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::max()const
{
  return applyUnaryReductionFunctor(PlayaFunctors::Max<Scalar>(this->comm()));
}


//===========================================================================
template <class Scalar> inline 
Scalar Vector<Scalar>::min()const
{
  return applyUnaryReductionFunctor(PlayaFunctors::Min<Scalar>(this->comm()));
}


//===========================================================================

template <class Scalar> inline 
const Scalar& Vector<Scalar>::operator[](int localIndex) const
{
  const SingleChunkVector<Scalar>* scv 
    = dynamic_cast<const SingleChunkVector<Scalar>* >(this->ptr().get());

  if (scv)
  {
    if (Debug::on)
    {
      PLAYA_BOUNDSCHECK(localIndex, 0, scv->chunkSize(), 
        "const Vector::operator[]()");
    }
    return scv->dataPtr()[localIndex];
  }
  
  int chunkBase = 0;
  while(this->ptr()->hasMoreChunks())
  {
    ConstDataChunk<Scalar> chunk = this->nextConstChunk();
    int chunkSize = chunk.size();
    if (localIndex >= chunkBase && localIndex < chunkBase+chunkSize)
    {
      this->ptr()->rewind();
      return chunk.values()[localIndex-chunkBase];
    }
    chunkBase += chunkSize;
  }
  
  TEST_FOR_EXCEPTION(true, std::runtime_error,
    "Vector operator[] local index " 
    << localIndex << " out of range [0, " << chunkBase << ")");

  return scv->dataPtr()[0]; // -Wall, will never be called.
}

//===========================================================================

template <class Scalar> inline 
Scalar& Vector<Scalar>::operator[](int localIndex)
{
  SingleChunkVector<Scalar>* scv 
    = dynamic_cast<SingleChunkVector<Scalar>* >(this->ptr().get());

  if (scv)
  {
    if (Debug::on)
    {
      PLAYA_BOUNDSCHECK(localIndex, 0, scv->chunkSize(), 
        "non-const Vector::operator[]()");
    }
    return scv->dataPtr()[localIndex];
  }
  
  int chunkBase = 0;
  while(this->ptr()->hasMoreChunks())
  {
    NonConstDataChunk<Scalar> chunk = this->nextChunk();
    int chunkSize = chunk.size();
    if (localIndex >= chunkBase && localIndex < chunkBase+chunkSize)
    {
      this->ptr()->rewind();
      return chunk.values()[localIndex-chunkBase];
    }
    chunkBase += chunkSize;
  }
  
  TEST_FOR_EXCEPTION(true, std::runtime_error,
    "Vector operator[] local index " 
    << localIndex << " out of range [0, " << chunkBase << ")");

  return scv->dataPtr()[0]; // -Wall, will never be called.
}

//===========================================================================

template <class Scalar> inline 
const Scalar& Vector<Scalar>::operator()(const BlockIterator<Scalar>& b,
  int localIndexWithinBlock) const
{
  return this->getBlock(b)[localIndexWithinBlock];
}


//===========================================================================

template <class Scalar> inline 
Scalar& Vector<Scalar>::operator()(const BlockIterator<Scalar>& b,
  int localIndexWithinBlock) 
{
  return this->getNonConstBlock(b)[localIndexWithinBlock];
}


//===========================================================================

template <class Scalar> inline
LoadableVector<Scalar>* loadable(Vector<Scalar> vec) 
{
  LoadableVector<Scalar>* lv 
    = dynamic_cast<LoadableVector<Scalar>* >(vec.ptr().get());
  TEST_FOR_EXCEPT(lv==0);
  return lv;
}

}




#endif
