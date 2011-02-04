/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_DEFAULTBLOCKVECTORSPACEIMPL_HPP
#define PLAYA_DEFAULTBLOCKVECTORSPACEIMPL_HPP

#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaDefaultBlockVectorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaExceptions.hpp"

namespace Playa
{

template <class Scalar> inline
DefaultBlockVectorSpace<Scalar>::
DefaultBlockVectorSpace(const Array<VectorSpace<Scalar> >& blocks) 
  : blocks_(blocks), baseGNI_(-1) 
{
  baseGNI_ = this->accumulateBaseGNI();
}

template <class Scalar> inline
const VectorSpace<Scalar>&  
DefaultBlockVectorSpace<Scalar>::getBlock(int b) const 
{
  TEST_FOR_EXCEPTION(b < 0 || b >= this->numBlocks(),
    RuntimeError, "block index b=" << b << " into vector space "
    << this->description() << " out of range [0,"
    << this->numBlocks() << ")");
  return blocks_[b];
}


template <class Scalar> inline
RCP<VectorBase<Scalar> > DefaultBlockVectorSpace<Scalar>
::createMember(const VectorSpace<Scalar>& self) const
{
  TEST_FOR_EXCEPTION(this != self.ptr().get(), RuntimeError,
    "inconsistent self-reference in DefaultBlockVectorSpace::"
    "createMember()");
  return rcp(new DefaultBlockVector<Scalar>(self));
}


template <class Scalar> inline
const MPIComm& DefaultBlockVectorSpace<Scalar>::comm() const
{
  return getBlock(0).comm();
}

template <class Scalar> inline
int DefaultBlockVectorSpace<Scalar>::baseGlobalNaturalIndex() const
{
  return baseGNI_;
}




/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2)
{
  Array<VectorSpace<Scalar> > x(2);
  x[0] = v1;
  x[1] = v2;
  return blockSpace<Scalar>(x);
}

/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2,
  const VectorSpace<Scalar>& v3)
{
  Array<VectorSpace<Scalar> > x(3);
  x[0] = v1;
  x[1] = v2;
  x[2] = v3;
  return blockSpace<Scalar>(x);
}

/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2,
  const VectorSpace<Scalar>& v3,
  const VectorSpace<Scalar>& v4)
{
  Array<VectorSpace<Scalar> > x(4);
  x[0] = v1;
  x[1] = v2;
  x[2] = v3;
  x[3] = v4;
  return blockSpace<Scalar>(x);
}

/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(const Array<VectorSpace<Scalar> >& x)
{
  RCP<const VectorSpaceBase<Scalar> > rtn 
    = rcp(new DefaultBlockVectorSpace<Scalar>(x));
  return rtn;
}


}


 

#endif
