/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_SERIAL_VECTOR_HPP
#define PLAYA_SERIAL_VECTOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaSingleChunkVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "Teuchos_Array.hpp"

namespace Playa
{
using namespace Teuchos;

template <class Scalar> class Vector;

/**
 * Playa implementation of a serial vector, implementing the LoadableVector
 * interface allowing an application to access elements. 
 * If created in SPMD, this will be replicated on
 * all processors.
 */
class SerialVector : public SingleChunkVector<double>,
                     public Describable
{
public:

  /** Construct with a vector space. */
  SerialVector(const VectorSpace<double>& vs);

  /** \name VectorBase interface */
  //@{
  /** Access to the space in which this vector lives */
  RCP<const VectorSpaceBase<double> > space() const {return vecSpace_.ptr();}
  //@}

  /** \name LoadableVector interface */
  //@{
  /** set a single element */
  void setElement(int globalIndex, const double& value);

  /** add to a single element */
  void addToElement(int globalIndex, const double& value);

  /** set a group of elements */
  void setElements(int numElems, const int* globalIndices, 
    const double* values);


  /** add to a group of elements */
  void addToElements(int numElems, const int* globalIndices, 
    const double* values);

  /** */
  void finalizeAssembly();
  //@}

  /** \name Diagnostics */
  //@{
  /** */
  std::string description() const ;
  //@}



  /** \name Access through global indices */
  //@{
  /** get the batch of elements at the given global indices */
  void getElements(const int* globalIndices, int numElems,
    Array<double>& elems) const ;
  //@}

  /** */
  static const SerialVector* getConcrete(const Vector<double>& x);
  /** */
  static SerialVector* getConcrete(Vector<double>& x);

      
 /** \name Single chunk data access interface */
  //@{
  /** */
  virtual const double* dataPtr() const {return &(data_[0]);}
  /** */
  virtual double* dataPtr() {return &(data_[0]);}

  /** Size of the (single) chunk of data values */
  virtual int chunkSize() const {return dim_;}
  //@}



  
private:

  VectorSpace<double> vecSpace_;

  Array<double> data_;

  int dim_;
};
  
}


#endif
