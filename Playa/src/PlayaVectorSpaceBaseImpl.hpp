/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTOR_SPACE_BASE_IMPL_HPP
#define PLAYA_VECTOR_SPACE_BASE_IMPL_HPP

#include "PlayaVectorSpaceBaseDecl.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_Array.hpp"

namespace Playa
{

template <class Scalar> inline
int VectorSpaceBase<Scalar>::accumulateBaseGNI() const 
{
  int nLocal = this->numLocalElements();
  Teuchos::Array<int> sums;
  int total = 0;
  Teuchos::MPIContainerComm<int>::accumulate(nLocal, sums, 
    total, this->comm());
  int rank = this->comm().getRank();
  return sums[rank];
}




}

#endif
