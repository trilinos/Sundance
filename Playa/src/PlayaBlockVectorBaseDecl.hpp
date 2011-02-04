/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_BLOCKVECTORBASEDECL_HPP
#define PLAYA_BLOCKVECTORBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorBaseDecl.hpp"
#include "Teuchos_Describable.hpp"

namespace Playa
{

template <class Scalar> class Vector;
using Teuchos::Describable;

/** 
 * Base class for blocked vectors 
 */
template <class Scalar>
class BlockVectorBase : public VectorBase<Scalar>,
                        public Describable
{
public:
  /** */
  BlockVectorBase() : currentBlock_(0) {}

  /** */
  virtual ~BlockVectorBase(){}

  /** */
  virtual void setBlock(int b, const Vector<Scalar>& block) = 0 ;

  /** */
  virtual const Vector<Scalar>& getBlock(int b) const = 0 ;

  /** */
  virtual Vector<Scalar> getNonConstBlock(int b) = 0 ;

  /** */
  virtual ConstDataChunk<Scalar> nextConstChunk() const ;

  /** */
  virtual NonConstDataChunk<Scalar> nextChunk() ;

  /** */
  virtual bool hasMoreChunks() const ;

  /** */
  virtual void rewind() const ;

  /** */
  virtual std::string description() const ;
  
  
private:
  mutable int currentBlock_;
};







/** */





      
      
  
  
}

#endif
