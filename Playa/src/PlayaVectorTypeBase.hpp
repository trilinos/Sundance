/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_VECTORTYPEBASE_HPP
#define PLAYA_VECTORTYPEBASE_HPP

#include "PlayaHandle.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp" 
#include "PlayaMatrixFactory.hpp" 
#include "PlayaGhostImporter.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 *
 */
template <class Scalar>
class VectorTypeBase
{
public:
  /** Virtual dtor */
  virtual ~VectorTypeBase() {;}

  /** create a distributed vector space.
   * @param dimension the dimension of the space 
   * @param nLocal number of indices owned by the local processor
   * @param locallyOwnedIndices array of indices owned by this processor  
   */
  virtual RCP<const VectorSpaceBase<Scalar> >
  createSpace(int dimension, 
    int nLocal,
    const int* locallyOwnedIndices,
    const MPIComm& comm) const = 0 ;
   

  /** Default implementation creates a vector space having 
   * nLocal elements on each processor. Serial types should override this
   * to produce a replicated space. */
  virtual VectorSpace<Scalar> 
  createEvenlyPartitionedSpace(const MPIComm& comm,
    int nLocal) const ;

  /**  
   * Create an importer for accessing ghost elements.
   * @param space the distributed vector space on which ghost elements
   * are to be shared
   * @param nGhost number of ghost elements needed by this processor
   * @param ghostIndices read-only C array of off-processor indices needed
   * by this processor.
   * @return A RCP to a GhostImporter object.
   */
  virtual RCP<GhostImporter<Scalar> > 
  createGhostImporter(const VectorSpace<Scalar>& space,
    int nGhost,
    const int* ghostIndices) const = 0 ;

    
  /**
   * Create a matrix factory of type compatible with this vector type,
   * sized according to the given domain and range spaces.
   */
  virtual RCP<MatrixFactory<Scalar> >
  createMatrixFactory(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range) const = 0 ;
    
    
};



/* Default implementation */
template <class Scalar> inline 
VectorSpace<Scalar> VectorTypeBase<Scalar>
::createEvenlyPartitionedSpace(const MPIComm& comm,
  int nLocal) const
{
  int rank = comm.getRank();
  int nProc = comm.getNProc();
  int dimension = nLocal * nProc;
  Array<int> locallyOwnedIndices(nLocal);
  int lowestLocalRow = rank*nLocal;
  for (int i=0; i<nLocal; i++)
  {
    locallyOwnedIndices[i] = lowestLocalRow + i;
  }
  return this->createSpace(dimension, nLocal, &(locallyOwnedIndices[0]), comm);
}

  
}

#endif
