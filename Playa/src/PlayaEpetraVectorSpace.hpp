/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_EPETRAVECTORSPACE_HPP
#define PLAYA_EPETRAVECTORSPACE_HPP

#include "PlayaDefs.hpp"
#include "Epetra_Map.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"

namespace Playa
{
using namespace Teuchos;


/**
 * Adaptor wrapping Epetra map in the Thyra vector space system.
 */
class EpetraVectorSpace : public VectorSpaceBase<double>
{
public:
  /** */
  EpetraVectorSpace(const RCP<const Epetra_Map>& map);

  /** */
  RCP<VectorBase<double> > 
  createMember(const VectorSpace<double>& self) const  ;
    
  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Public overridden from VectorSpace */
  //@{

  /** */
  int dim() const {return globalDim_;}
  
  /** */
  int numLocalElements() const {return numLocalElements_;}
  
  /** */
  int baseGlobalNaturalIndex() const {return baseGlobalNaturalIndex_;}

  /** */
  bool isCompatible(const VectorSpaceBase<double>* other) const ;

  int numBlocks() const {return 1;}

  //@}

 /** */
  const RCP<const Epetra_Map>& epetraMap() const 
    {return epetraMap_;}


  /** */
  virtual const MPIComm& comm() const {return comm_;}

protected:

  MPIComm epetraCommToTeuchosMPIComm(const Epetra_Comm& epComm) ;
  
private:
  RCP<const Epetra_Map> epetraMap_;

  MPIComm comm_;

  int globalDim_;

  int baseGlobalNaturalIndex_;

  int numLocalElements_;
};
  
}

#endif
