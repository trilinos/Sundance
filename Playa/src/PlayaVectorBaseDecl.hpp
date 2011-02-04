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

};


  
  
}

#endif
