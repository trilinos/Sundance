/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;


MatrixLaplacian1D::MatrixLaplacian1D(int nLocalRows, 
  const VectorType<double>& type)
  : OperatorBuilder<double>(nLocalRows, type), op_()
{
  init(nLocalRows, type);
}




void MatrixLaplacian1D::init(int nLocalRows, 
  const VectorType<double>& type)
{
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  RCP<MatrixFactory<double> > mFact;
  mFact = vecType().createMatrixFactory(domain(), range());

  if (domain().dim() == domain().numLocalElements())
  {
    rank = 0;
    nProc = 1;
  }

  int lowestLocalRow = nLocalRows * rank;

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(mFact.get());
  if (icmf)
  {
    for (int i=0; i<nLocalRows; i++)
    {
      int row = lowestLocalRow + i;
      Array<int> colIndices;
      if (rank==0 && i==0)
      {
        colIndices = tuple(row, row+1);
      }
      else if (rank==nProc-1 && i==nLocalRows-1)
      {
        colIndices = tuple(row-1, row);
      }
      else
      {
        colIndices = tuple(row-1, row, row+1);
      }
      icmf->initializeNonzerosInRow(row, colIndices.size(),
        &(colIndices[0]));
    }
    icmf->finalize();
  }
      
  op_ = mFact->createMatrix();

  double h = 1.0/((double) domain().dim() + 1);
      
  RCP<LoadableMatrix<double> > mat = op_.matrix();

  /* fill in with the Laplacian operator */
  for (int i=0; i<nLocalRows; i++)
  {
    int row = lowestLocalRow + i;
    Array<int> colIndices;
    Array<double> colVals;
    if (rank==0 && i==0)
    {
      colIndices = tuple(row, row+1);
      colVals = tuple(2.0/h/h, -1.0/h/h);
    }
    else if (rank==nProc-1 && i==nLocalRows-1)
    {
      colIndices = tuple(row-1, row);
      colVals = tuple(-1.0/h/h, 2.0/h/h);
    }
    else
    {
      colIndices = tuple(row-1, row, row+1);
      colVals = tuple(-1.0/h/h, 2.0/h/h, -1.0/h/h);
    }
    mat->addToRow(row, colIndices.size(), 
      &(colIndices[0]), &(colVals[0]));
  }
}
