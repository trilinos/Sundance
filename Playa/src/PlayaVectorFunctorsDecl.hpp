/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORFUNCTORDECL_HPP
#define PLAYA_VECTORFUNCTORDECL_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_RCP.hpp"
#include "PlayaGeneralizedIndex.hpp"

namespace Playa
{

using Teuchos::RCP;
using Teuchos::MPIComm;

/**
 * This traits class specifies the return type of a reduction functor. 
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar, class FunctorType>
class VectorFunctorTraits
{
public:
  typedef Scalar ReturnType;
};



/** */
template <class Scalar>
class ReductionFunctorBase
{
public:
  /** */
  ReductionFunctorBase(
    const MPIComm& comm
    )
    : comm_(comm), gi_(rcp(new GeneralizedIndex())) {}

  /** */
  virtual void postProc() const = 0 ;

  /** */
  const MPIComm& comm() const {return comm_;}

  /** */
  RCP<GeneralizedIndex> currentGI() const {return gi_;}

private:
  MPIComm comm_;
  mutable RCP<GeneralizedIndex> gi_;
};

  
  
}

#endif
