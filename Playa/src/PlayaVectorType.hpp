/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_VECTORTYPE_HPP
#define PLAYA_VECTORTYPE_HPP

#include "PlayaHandle.hpp"
#include "PlayaVectorTypeBase.hpp"
#include "PlayaVectorSpaceDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * Vector type objects are used by the application code to create
 * vector spaces and operators of a given type.
 */
template <class Scalar>
class VectorType : public Playa::Handle<VectorTypeBase<Scalar> >
{
public:
  HANDLE_CTORS(VectorType<Scalar>, VectorTypeBase<Scalar>);
   

  /** create a vector space having nLocal elements on each processor */
  VectorSpace<Scalar> createEvenlyPartitionedSpace(const MPIComm& comm,
    int nLocal) const ;

  /** create a distributed vector space.
   * @param dimension the dimension of the space
   * @param nLocal number of indices owned by the local processor
   * @param locallyOwnedIndices array of indices owned by this processor  
   */
  VectorSpace<Scalar> createSpace(int dimension, 
    int nLocal,
    const int* locallyOwnedIndices,
    const MPIComm& comm) const ;



  /** 
   * Create an importer for ghost elements
   **/
  RCP<GhostImporter<Scalar> > 
  createGhostImporter(const VectorSpace<Scalar>& space,
    int nGhost,
    const int* ghostIndices) const ;

  /**
   * Create a matrix factory of type compatible with this vector type,
   * sized according to the given domain and range spaces.
   */
  virtual RCP<MatrixFactory<Scalar> >
  createMatrixFactory(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range) const ;
                                                      
};




template <class Scalar> inline 
VectorSpace<Scalar> VectorType<Scalar>::createSpace(int dimension,
  int nLocal,
  const int* locallyOwnedIndices, const MPIComm& comm) const
{
  return this->ptr()->createSpace(dimension, nLocal, locallyOwnedIndices, comm);
}

template <class Scalar> inline 
VectorSpace<Scalar> VectorType<Scalar>
::createEvenlyPartitionedSpace(const MPIComm& comm,
  int nLocal) const
{
  return this->ptr()->createEvenlyPartitionedSpace(comm, nLocal);
}


template <class Scalar> inline 
RCP<GhostImporter<Scalar> > 
VectorType<Scalar>::createGhostImporter(const VectorSpace<Scalar>& space,
  int nGhost,
  const int* ghostIndices) const
{
  return this->ptr()->createGhostImporter(space, nGhost, ghostIndices);
}

template <class Scalar> inline
RCP<MatrixFactory<Scalar> >
VectorType<Scalar>::createMatrixFactory(const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range) const
{
  return this->ptr()->createMatrixFactory(domain, range);
}

  
}

#endif
