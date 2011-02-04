/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_SERIAL_VECTORSPACE_HPP
#define PLAYA_SERIAL_VECTORSPACE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * 
 */
class SerialVectorSpace : public VectorSpaceBase<double>
{
public:

  /** */
  SerialVectorSpace(int dim);
    

  /** @name Overridden from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** */
  RCP<VectorBase<double> > createMember(const VectorSpace<double>& self) const ;
  
  /** */
  int dim() const {return dim_;}

  /** */
  int numBlocks() const {return 1;}
  
  /** */
  int numLocalElements() const {return dim();}
  
  /** */
  int baseGlobalNaturalIndex() const {return 0;}
  
  /** */
  bool isCompatible(const VectorSpaceBase<double>* other) const ;

  /** */
  virtual const MPIComm& comm() const {return comm_;}

private:
  int dim_;
  MPIComm comm_;
};
  
}

#endif
