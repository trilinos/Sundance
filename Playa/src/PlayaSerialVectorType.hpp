/* @HEADER@ */
//
/* @HEADER@ */

#ifndef PLAYA_SERIAL_VECTORTYPE_HPP
#define PLAYA_SERIAL_VECTORTYPE_HPP

#include "PlayaSerialVectorSpace.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaVectorTypeBase.hpp"
//#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;
  
/**
 * \!brief Serial vector type is a factory for serial vector spaces. If used
in a SPMD program the spece will be replicated on each processor. 
 */
class SerialVectorType : public VectorTypeBase<double>,
                         public Handleable<VectorTypeBase<double> >,
                         public Printable,
                         public Describable
{
public:
  /** Construct a vector type */
  SerialVectorType();
      
  /** virtual dtor */
  virtual ~SerialVectorType() {;}

  /** create a distributed vector space.
   * @param dimension the dimension of the space 
   * @param nLocal number of indices owned by the local processor
   * @param locallyOwnedIndices array of indices owned by this processor  
   */
  RCP<const VectorSpaceBase<double> > 
  createSpace(int dimension, 
    int nLocal,
    const int* locallyOwnedIndices,
    const MPIComm& comm) const ;


  /** Default implementation creates a vector space having 
   * nLocal elements on each processor. Serial types should override this
   * to produce a replicated space. */
  virtual VectorSpace<double> 
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
  RCP<GhostImporter<double> > 
  createGhostImporter(const VectorSpace<double>& space,
    int nGhost,
    const int* ghostIndices) const ;


  /**
   * Create a matrix factory of type compatible with this vector type,
   * sized according to the given domain and range spaces.
   */
  RCP<MatrixFactory<double> >
  createMatrixFactory(const VectorSpace<double>& domain,
    const VectorSpace<double>& range) const ;

    

  /** \name Printable interface */
  //@{
  /** Print to stream */
  void print(std::ostream& os) const {os << description();}
  //@}

  GET_RCP(VectorTypeBase<double>);
};
  
}

#endif
