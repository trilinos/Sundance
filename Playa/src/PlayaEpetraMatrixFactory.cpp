/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaEpetraMatrixFactory.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  
#include "PlayaVectorDecl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaLinearOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaVectorImpl.hpp"
#endif


namespace Playa
{

using namespace Teuchos;


EpetraMatrixFactory::EpetraMatrixFactory(
  const RCP<const EpetraVectorSpace>& domain,
  const RCP<const EpetraVectorSpace>& range)
  : graph_(rcp(new Epetra_CrsGraph(Copy, *(range->epetraMap()), 0))),
    range_(range),
    domain_(domain)
{}


void EpetraMatrixFactory::finalize()
{
  int ierr = graph_->FillComplete(*(domain_->epetraMap()), *(range_->epetraMap()));

  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
    "EpetraMatrixFactory::finalize() failed during call "
    "to FillComplete(). Error code was " << ierr);

  if (!graph_->StorageOptimized())
  {
    ierr = graph_->OptimizeStorage();
      
    TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
      "EpetraMatrixFactory::freezeValues() failed during call "
      "to OptimizeStorage(). Error code was " << ierr);
  }
}

void EpetraMatrixFactory::initializeNonzerosInRow(int globalRowIndex,
  int nElemsToInsert,
  const int* globalColumnIndices)
{
  int ierr = graph_->InsertGlobalIndices(globalRowIndex,
    nElemsToInsert,
    (int*) globalColumnIndices);
  
  TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
    "failed to add to row " << globalRowIndex
    << " in EpetraMatrixFactory::setRowValues() with nnz="
    << nElemsToInsert 
    << ". Error code was " << ierr);
}


void EpetraMatrixFactory::initializeNonzeroBatch(int numRows, 
  int rowBlockSize,
  const int* globalRowIndices,
  int numColumnsPerRow,
  const int* globalColumnIndices,
  const int* skipRow)
{
  int numRowBlocks = numRows/rowBlockSize;
  int row = 0;
  
  for (int rb=0; rb<numRowBlocks; rb++)
  {
    const int* cols = globalColumnIndices + rb*numColumnsPerRow;
    for (int r=0; r<rowBlockSize; r++, row++)
    {
      if (skipRow[row]) continue;
      graph_->InsertGlobalIndices(globalRowIndices[row], 
        numColumnsPerRow, (int*) cols);
    }
  }
}


void EpetraMatrixFactory::configure(int lowestRow,
  const std::vector<int>& rowPtrs,
  const std::vector<int>& nnzPerRow,
  const std::vector<int>& data)
{

  graph_ = rcp(new Epetra_CrsGraph(Copy, *(range_->epetraMap()),
      (const int*) &(nnzPerRow[0]),
      true));
  
  for (unsigned int i=0; i<rowPtrs.size(); i++)
  {
    graph_->InsertGlobalIndices(lowestRow + i, nnzPerRow[i],
      (int*) &(data[rowPtrs[i]]));
  }

  finalize();
}

const Epetra_CrsGraph& EpetraMatrixFactory::graph() const 
{
  return *(graph_.get());
}


LinearOperator<double> EpetraMatrixFactory::createMatrix() const
{
  RCP<const VectorSpaceBase<double> > dp = epDomain();
  RCP<const VectorSpaceBase<double> > rp = epRange();
  VectorSpace<double> d = dp;
  VectorSpace<double> r = rp;

  RCP<LinearOperatorBase<double> > A 
    = rcp(new EpetraMatrix(graph(), d, r));
  return A;
}

}
