/* @HEADER@ */
//
 /* @HEADER@ */

#include "PlayaVectorDecl.hpp"
#include "PlayaSerialVector.hpp"
#include "PlayaSerialVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Teuchos;
using namespace Playa;

SerialVector::SerialVector(const VectorSpace<double>& vs)
  : SingleChunkVector<double>(),
    vecSpace_(vs),
    data_(vs.dim()),
    dim_(vs.dim())
{
  const SerialVectorSpace* rvs 
    = dynamic_cast<const SerialVectorSpace*>(vs.ptr().get());
  TEST_FOR_EXCEPTION(rvs==0, std::runtime_error,
    "could not cast vector space to SerialVectorSpace in "
    "SerialVector ctor");
}


void SerialVector::setElement(OrdType index, const double& value)
{
  data_[index] = value;
}

void SerialVector::addToElement(OrdType index, const double& value)
{
  data_[index] += value;
}

void SerialVector::setElements(int numElems, const int* globalIndices,
  const double* values)
{
  for (int i=0; i<numElems; i++)
  {
    data_[globalIndices[i]] = values[i];
  }
}

void SerialVector::addToElements(int numElems, const int* globalIndices,
  const double* values)
{
  for (int i=0; i<numElems; i++)
  {
    data_[globalIndices[i]] += values[i];
  }
}

const SerialVector* SerialVector::getConcrete(const Vector<double>& x)
{
  const SerialVector* rtn = dynamic_cast<const SerialVector*>(x.ptr().get());
  TEST_FOR_EXCEPT(rtn==0);
  return rtn;
}

SerialVector* SerialVector::getConcrete(Vector<double>& x)
{
  SerialVector* rtn = dynamic_cast<SerialVector*>(x.ptr().get());
  TEST_FOR_EXCEPT(rtn==0);
  return rtn;
}

void SerialVector::finalizeAssembly()
{
  // no-op
}

void SerialVector::getElements(const OrdType* globalIndices, OrdType numElems,
  Array<double>& elems) const
{
  elems.resize(numElems);
  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*this)[globalIndices[i]];
  }
}

std::string SerialVector::description() const 
{
  std::ostringstream oss;
  oss << "SerialVector[dim=" << dim_ << "]" ;
  return oss.str();
}

