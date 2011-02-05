/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SINGLECHUNKVECTOR_HPP
#define PLAYA_SINGLECHUNKVECTOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorBaseDecl.hpp"

namespace Playa
{
/**
 * Base class for vector types that have all on-processor data in a single
 * contiguous chunk
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class SingleChunkVector : public VectorBase<Scalar>
{
public:
  /** */
  SingleChunkVector() : rewound_(true) {}
  /** virtual dtor */
  virtual ~SingleChunkVector() {}

  /** */
  virtual ConstDataChunk<Scalar> nextConstChunk() const 
    {rewound_ = false; return ConstDataChunk<Scalar>(chunkSize(), dataPtr());}

  /** */
  virtual NonConstDataChunk<Scalar> nextChunk() 
    {rewound_ = false; return NonConstDataChunk<Scalar>(chunkSize(), dataPtr());}

  /** */
  virtual bool hasMoreChunks() const 
    {return rewound_;}

  /** */
  virtual void rewind() const 
    {rewound_=true;}

  /** \name Access to local elements */
  //@{
  /** read the element at the given local index */
  virtual const double& operator[](int localIndex) const 
    {return dataPtr()[localIndex];}

  /** writable access to the element at the given local index */
  virtual double& operator[](int localIndex) 
    {return dataPtr()[localIndex];}
  //@}

  /** */
  virtual int chunkSize() const = 0 ;
  /** */
  virtual const Scalar* dataPtr() const = 0 ;
  /** */
  virtual Scalar* dataPtr() = 0 ;

private:
  mutable bool rewound_;
};
}


#endif
