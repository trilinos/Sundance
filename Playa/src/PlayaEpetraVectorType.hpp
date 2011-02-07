/* @HEADER@ */
//
/* @HEADER@ */

#ifndef PLAYA_EPETRAVECTORTYPE_HPP
#define PLAYA_EPETRAVECTORTYPE_HPP

#include "PlayaEpetraVectorSpace.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaVectorTypeBase.hpp"


namespace Playa
{
using namespace Teuchos;
  
/**
 * \!brief Epetra vector type is a factory for epetra vector spaces
 */
class EpetraVectorType : public VectorTypeBase<double>,
                         public Playa::Handleable<VectorTypeBase<double> >,
                         public Printable,
                         public Describable
{
public:
  /** Construct a vector type */
  EpetraVectorType();
      
  /** virtual dtor */
  virtual ~EpetraVectorType() {;}

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
