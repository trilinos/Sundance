/* @HEADER@ */
//
 /* @HEADER@ */

#include "PlayaSerialVectorSpace.hpp"
#include "PlayaSerialVector.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaOut.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorSpaceImpl.hpp"
#endif




namespace Playa
{
using Teuchos::RCP;
using Teuchos::MPIComm;

SerialVectorSpace::SerialVectorSpace(int dim)
  : dim_(dim), comm_(MPIComm::self())
{}


// Overridden from VectorSpaceBase

RCP<VectorBase<double> >
SerialVectorSpace::createMember(const VectorSpace<double>& self) const
{
  TEST_FOR_EXCEPTION(self.ptr().get() != this,
    InternalError, 
    "inconsistency between space and self-reference in SerialVectorSpace::createMember()");

  return rcp(new SerialVector(self));
}

bool SerialVectorSpace::isCompatible(const VectorSpaceBase<double>* other) const
{
  const SerialVectorSpace* svs = dynamic_cast<const SerialVectorSpace*>(other);
  if (svs == 0) return false;
  return this->dim() == svs->dim();
}

string SerialVectorSpace::description() const
{
  std::string rtn = "SerialVS[d=" + Teuchos::toString(this->dim()) + "]";
  return rtn;
}

}


