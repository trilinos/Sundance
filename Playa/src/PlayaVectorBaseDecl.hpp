/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORBASEDECL_HPP
#define PLAYA_VECTORBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"
#include "PlayaDataChunk.hpp"

namespace Playa
{
template <class Scalar> class Vector;

/**
 * 
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class VectorBase
{
public:
  /** virtual dtor */
  virtual ~VectorBase() {;}

  /** Access to the space in which this vector lives */
  virtual RCP<const VectorSpaceBase<Scalar> > space() const = 0 ;

  /** */
  virtual ConstDataChunk<Scalar> nextConstChunk() const = 0 ;

  /** */
  virtual NonConstDataChunk<Scalar> nextChunk() = 0 ;

  /** */
  virtual bool hasMoreChunks() const = 0 ;

  /** */
  virtual void rewind() const = 0 ;

  /** */
  virtual int numBlocks() const {return 1;}

  /** Carry out the operation
   * (*this) = gamma*(*this) + alpha*x ;
   */
  virtual void update(const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& gamma) = 0;

  /** Carry out the operation
   * (*this) = gamma*(*this) + alpha*x + beta*y;
   */
  virtual void update(
    const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& beta, const VectorBase<Scalar>* y,
    const Scalar& gamma) = 0 ;

  /** */
  virtual void update(
    const Scalar& alpha, const VectorBase<Scalar>* x,
    const Scalar& beta, const VectorBase<Scalar>* y,
    const Scalar& gamma, const VectorBase<Scalar>* z,
    const Scalar& delta) = 0 ;

  /** */
  virtual Scalar dot(const VectorBase<Scalar>* other) const = 0 ;

  /** */
  virtual Scalar norm2() const = 0 ;
    
};


  
  
}

#endif
