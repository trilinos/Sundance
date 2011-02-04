/* @HEADER@ */
//
 /* @HEADER@ */

#include "PlayaEpetraVectorSpace.hpp"
#include "PlayaEpetraVector.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "PlayaOut.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "Teuchos_DefaultMpiComm.hpp"
#endif


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorSpaceImpl.hpp"
#endif




using Teuchos::RCP;

namespace Playa
{

EpetraVectorSpace::EpetraVectorSpace(const RCP<const Epetra_Map>& m)
  : epetraMap_(m),
    comm_(epetraCommToTeuchosMPIComm(m->Comm())),
    globalDim_(m->NumGlobalElements()),
    baseGlobalNaturalIndex_(m->MinMyGID()),
    numLocalElements_(m->NumMyElements())
{}

MPIComm EpetraVectorSpace::epetraCommToTeuchosMPIComm(
  const Epetra_Comm& epComm) 
{
#ifdef HAVE_MPI
  const Epetra_MpiComm* c1 = dynamic_cast<const Epetra_MpiComm*>(&epComm);
  if (c1) return MPIComm(c1->Comm());
#endif
  return MPIComm::self();
}

bool EpetraVectorSpace::isCompatible(const VectorSpaceBase<double>* other) const 
{
  const EpetraVectorSpace* evs = dynamic_cast<const EpetraVectorSpace*>(other);
  if (evs == 0) return false;
  return epetraMap_->SameAs(*(evs->epetraMap_));
}

// Overridden from VectorSpace

Teuchos::RCP<VectorBase<double> >
EpetraVectorSpace::createMember(const VectorSpace<double>& self) const
{
  TEST_FOR_EXCEPTION(self.ptr().get() != this,
    InternalError, 
    "inconsistency between space and self-reference in EpetraVectorSpace::createMember()");
  return rcp(new EpetraVector(self));
}

string EpetraVectorSpace::description() const
{
  std::string rtn = "EpetraVS[d=" + Teuchos::toString(dim());
  if (numLocalElements() != dim()) rtn += ", local="
    + Teuchos::toString(numLocalElements());
  rtn += "]";
  return rtn;
}

}



