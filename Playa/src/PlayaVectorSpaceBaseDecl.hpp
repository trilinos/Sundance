/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_VECTORSPACEBASEDECL_HPP
#define PLAYA_VECTORSPACEBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Describable.hpp"

namespace Teuchos
{
class MPIComm;
}

namespace Playa
{
using Teuchos::RCP;
using Teuchos::MPIComm;
using Teuchos::Describable;

template <class Scalar> class VectorSpace;
template <class Scalar> class VectorBase;
template <class Scalar> class BlockIterator;


/**
 * 
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar>
class VectorSpaceBase : public Describable
{
public:
  /** virtual dtor */
  virtual ~VectorSpaceBase() {;}

  /** */
  virtual RCP<VectorBase<Scalar> > createMember(const VectorSpace<Scalar>& self) const = 0 ;

  /** */
  virtual int dim() const = 0 ;

  /** */
  virtual int numLocalElements() const = 0 ;

  /** */
  virtual int baseGlobalNaturalIndex() const = 0 ;

  /** */
  virtual bool isCompatible(const VectorSpaceBase<Scalar>* other) const = 0 ;

  /** */
  virtual const MPIComm& comm() const = 0 ;

  /** */
  virtual int numBlocks() const = 0 ;

protected:
  /** */
  int accumulateBaseGNI() const ;
};

       
  
}

#endif
