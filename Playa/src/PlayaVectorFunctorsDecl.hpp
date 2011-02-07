/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORFUNCTORDECL_HPP
#define PLAYA_VECTORFUNCTORDECL_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_RCP.hpp"
#include "PlayaGeneralizedIndex.hpp"

namespace PlayaFunctors
{

using Playa::GeneralizedIndex;
using Teuchos::RCP;
using Teuchos::MPIComm;

/**
 * \brief This traits class specifies the return type of a reduction functor. 
 * If not specialized, the default return type will be a Scalar.
 *
 * @author Kevin Long (kevin.long@ttu.edu)
 */
template <class Scalar, class FunctorType>
class VectorFunctorTraits
{
public:
  typedef Scalar ReturnType;
};



/** 
 * \brief Base class for reduction functors
 *
 * @author Kevin Long (kevin.long@ttu.edu)
*/
template <class Scalar>
class ReductionFunctorBase
{
public:
  /** Construct with a communicator */
  ReductionFunctorBase(
    const MPIComm& comm
    )
    : comm_(comm) {}

  /** Callback for any postprocessing step (for example, MPI all-reduce) */
  virtual void postProc() const = 0 ;

protected:

  /** Return the MPI communicator */
  const MPIComm& comm() const {return comm_;}

private:
  MPIComm comm_;
};

  
/**
 * \brief IndexedValue is the return type for reduction operations such
 * as MinLoc that return a location and a value. 
 */
template <class Scalar>
struct IndexedValue
{
  /** Value */
  Scalar what;
  /** Index */
  int where;
};
  
}

#endif
