/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_DEFAULT_BLOCK_VECTOR_SPACE_DECL_HPP
#define PLAYA_DEFAULT_BLOCK_VECTOR_SPACE_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaBlockVectorSpaceDecl.hpp"
#include "Teuchos_Array.hpp"



namespace Playa
{
using Teuchos::Array;

/**
 * This is the default implementation of a blocked vector space, which
 * simply stores blocks in an array.
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class DefaultBlockVectorSpace : public BlockVectorSpaceBase<Scalar>
{
public:
  /** ctor */
  DefaultBlockVectorSpace(const Array<VectorSpace<Scalar> >& blocks) ;

  /** Get a block specified by an integer index. This function
   * should hrow an exception if the index is out of range */
  virtual const VectorSpace<Scalar>& getBlock(int b) const ;

  /** Create a block vector */
  virtual RCP<VectorBase<Scalar> > createMember(const VectorSpace<Scalar>& self) const ;

  /** Return the communicator */
  virtual const MPIComm& comm() const ;

  /** Return the lowest index of the first block */
  virtual int baseGlobalNaturalIndex() const ;

  /** */
  int numBlocks() const {return blocks_.size();}


private:
  Array<VectorSpace<Scalar> > blocks_;
  int baseGNI_;
};


/** \relates VectorSpace */
template <class Scalar> 
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1);

/** \relates VectorSpace */
template <class Scalar> 
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2);

/** \relates VectorSpace */
template <class Scalar> 
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2,
  const VectorSpace<Scalar>& v3);

/** \relates VectorSpace */
template <class Scalar> 
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2,
  const VectorSpace<Scalar>& v3,
  const VectorSpace<Scalar>& v4);

/** \relates VectorSpace */
template <class Scalar> 
VectorSpace<Scalar> blockSpace(const Array<VectorSpace<Scalar> >& x);

}

#endif
