/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LOADABLEVECTOR_HPP
#define PLAYA_LOADABLEVECTOR_HPP

#include "PlayaDefs.hpp"

namespace Playa
{
/**
 * LoadableVector defines an interface through which elements can 
 * be loaded into a vector. Element loading is used extensively
 * by application codes in creating vectors, 
 * but should never be used by high-performance solver codes.
 *
 * A vector type that will be
 * used in a context where loading is required should multiply inherit
 * from both VectorBase andLoadableVector.
 * 
 * Elements can by loaded one at a time
 * or in batches. The methods to load single elements arew pure virtual
 * and thus must be defined by derived classes. 
 * Loading in batches will usually be more efficient
 * provided the underlying vector implementation supports it. 
 * For those types not supporting batch loading, LoadableVector provides
 * default batch loading functions which delegate to single-element loading.
 *
 * Elements can by loaded either by setting a value, or adding to an 
 * existing value. The latter will typically by used in finite-element
 * codes.
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class LoadableVector 
{
public:
  /** virtual dtor */
  virtual ~LoadableVector() {;}

  /** set a single element at the given global index */
  virtual void setElement(int globalIndex, const Scalar& value) = 0 ;

  /** add to the existing value of 
   * a single element at the given global index */
  virtual void addToElement(int globalIndex, const Scalar& value) = 0 ;

  /** set a group of elements */
  virtual void setElements(int numElems, 
    const int* globalIndices, 
    const Scalar* values) ;

  /** add to a group of elements */
  virtual void addToElements(int numElems, 
    const int* globalIndices, 
    const Scalar* values);

  /** Do whatever finalization steps are needed by the implementation,
      for instance, synchronizing border elements. The default implementation
      * is a no-op. */
  virtual void finalizeAssembly() {;}
};

/* Default implementation of setElements makes multiple calls to
 * setElement(). If at all possible, this should be overridden
 * with a method specialized to the underlying type.  */
template <class Scalar> 
inline void LoadableVector<Scalar>::setElements(int numElems, 
  const int* globalIndices, 
  const Scalar* values)
{
  for (int i=0; i<numElems; i++)
  {
    setElement(globalIndices[i], values[i]);
  }
}

/* Default implementation of addToElements makes multiple calls to
 * addToElement(). If at all possible, this should be overridden
 * with a method specialized to the underlying type.  */
template <class Scalar> 
inline void LoadableVector<Scalar>::addToElements(int numElems, 
  const int* globalIndices, 
  const Scalar* values)
{
  for (int i=0; i<numElems; i++)
  {
    addToElement(globalIndices[i], values[i]);
  }
}

  
  
}

#endif
