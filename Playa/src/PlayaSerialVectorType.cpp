/* @HEADER@ */
//
 /* @HEADER@ */


#include "PlayaSerialVectorType.hpp"
#include "PlayaSerialVectorSpace.hpp"
#include "PlayaDenseSerialMatrixFactory.hpp"
#include "PlayaSerialGhostImporter.hpp"
#include "PlayaOut.hpp"

#include "Teuchos_RefCountPtr.hpp"

namespace Playa
{

using namespace Teuchos;

SerialVectorType::SerialVectorType()
{;}


RCP<const VectorSpaceBase<double> > 
SerialVectorType::createSpace(int dimension,
  int nLocal,
  const int* localIndices,
  const MPIComm& comm) const
{
  TEST_FOR_EXCEPTION(nLocal < 0, std::runtime_error, "negative vector size n=" << nLocal);
  TEST_FOR_EXCEPTION(dimension != nLocal, std::runtime_error, 
    "nLocal=" << nLocal << " and dimension=" << dimension
    << " should be equal for a replicated space");

	return rcp(new SerialVectorSpace(dimension));
}

RCP<GhostImporter<double> > 
SerialVectorType::createGhostImporter(const VectorSpace<double>& space,
                                      int nGhost,
                                      const int* ghostIndices) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const SerialVectorSpace*>(space.ptr().get())==0, std::runtime_error, "expected " 
    << space << " to be a SerialVectorSpace");
  return rcp(new SerialGhostImporter());
}

RCP<MatrixFactory<double> >
SerialVectorType::createMatrixFactory(const VectorSpace<double>& domain,
                                      const VectorSpace<double>& range) const
{
  RCP<MatrixFactory<double> > rtn 
    = rcp(new DenseSerialMatrixFactory(domain, range));

  return rtn;
}


VectorSpace<double> SerialVectorType
::createEvenlyPartitionedSpace(const MPIComm& /* comm */,
  int nLocal) const
{
  RCP<const VectorSpaceBase<double> > rtn = rcp(new SerialVectorSpace(nLocal));
  return rtn;
}

}




