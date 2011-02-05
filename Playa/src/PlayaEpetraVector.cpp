/* @HEADER@ */
//
 /* @HEADER@ */

#include "PlayaEpetraVector.hpp"
#include "PlayaEpetraVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using Teuchos::RCP;
using namespace Playa;

EpetraVector::EpetraVector(const VectorSpace<double>& vs)
  : SingleChunkVector<double>(), 
    vecSpace_(vs),
    epetraVec_(), 
    numLocalElements_(vs.numLocalElements())
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.ptr().get());
  TEST_FOR_EXCEPTION(epvs==0, std::runtime_error,
    "could not cast vector space to EpetraVectorSpace in "
    "EpetraVector ctor");

  epetraVec_ = rcp(new Epetra_Vector(*(epvs->epetraMap())));
}



EpetraVector
::EpetraVector(const VectorSpace<double>& vs,
  const RCP<Epetra_Vector>& vec)
  : SingleChunkVector<double>(), 
    vecSpace_(vs), 
    epetraVec_(vec), 
    numLocalElements_(vs.numLocalElements())
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.ptr().get());
  TEST_FOR_EXCEPTION(epvs==0, std::runtime_error,
    "could not cast vector space to EpetraVectorSpace in "
    "EpetraVector ctor");
}





double& EpetraVector::operator[](int globalIndex) 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(globalIndex)];
}

void EpetraVector::setElement(int index, const double& value)
{
  int loc_index[1] = { index };
  epetraVec()->ReplaceGlobalValues(1, const_cast<double*>(&value), 
    loc_index);
}

void EpetraVector::addToElement(int index, const double& value)
{
//  cout << "adding (" << index << ", " << value << ")" << std::endl;
  int loc_index[1] = { index };
  epetraVec()->SumIntoGlobalValues(1, const_cast<double*>(&value), 
    loc_index);
}

const double& EpetraVector::getElement(int index) const 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(index)];
}

void EpetraVector::getElements(const int* globalIndices, int numElems,
  Teuchos::Array<double>& elems) const
{
  elems.resize(numElems);
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  RCP<const Epetra_Vector> epv = epetraVec();

  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*epv)[myMap.LID(globalIndices[i])];
  }
}

void EpetraVector::setElements(int numElems, const int* globalIndices,
  const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->ReplaceGlobalValues(numElems, globalIndices, values);
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "ReplaceGlobalValues returned "
    "ierr=" << ierr << " in EpetraVector::setElements()");
}

void EpetraVector::addToElements(int numElems, const int* globalIndices,
  const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->SumIntoGlobalValues(numElems, globalIndices, values);
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "SumIntoGlobalValues returned "
    "ierr=" << ierr << " in EpetraVector::addToElements()");
}

void EpetraVector::finalizeAssembly()
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  vec->GlobalAssemble();
}


void EpetraVector::print(std::ostream& os) const 
{
  epetraVec()->Print(os);
}


const Epetra_Vector& EpetraVector::getConcrete(const Playa::Vector<double>& tsfVec)
{
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}

Epetra_Vector& EpetraVector::getConcrete(Playa::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}


Epetra_Vector* EpetraVector::getConcretePtr(Playa::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return epv->epetraVec().get();
}




