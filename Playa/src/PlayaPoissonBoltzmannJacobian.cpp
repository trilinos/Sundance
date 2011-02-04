/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaPoissonBoltzmannJacobian.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"
#include "PlayaTabs.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif
using namespace Playa;
using namespace Teuchos;


PoissonBoltzmannJacobian
::PoissonBoltzmannJacobian(int nLocalRows, 
  const VectorType<double>& type)
  : OperatorBuilder<double>(nLocalRows, type), op_(), nLocalRows_(nLocalRows),
    h_(1.0)
{
  h_ = 1.0/((double) domain().dim() - 1);
}

void PoissonBoltzmannJacobian::setEvalPoint(const Vector<double>& x)
{
  Tabs tab;
  Out::os() << tab << "in PBJ::setEvalPoint()" << std::endl;
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  RCP<MatrixFactory<double> > mFact 
    = vecType().createMatrixFactory(domain(), range());
  
  int lowestLocalRow = nLocalRows_ * rank;

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(mFact.get());
  for (int i=0; i<nLocalRows_; i++)
  {
    int row = lowestLocalRow + i;
    Array<int> colIndices;
    if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows_-1))
    {
      colIndices = tuple(row);
    }
    else
    {
      colIndices = tuple(row-1, row, row+1);
    }
    icmf->initializeNonzerosInRow(row, colIndices.size(),
      &(colIndices[0]));
  }
  icmf->finalize();
      
  op_ = mFact->createMatrix();
      
  RCP<LoadableMatrix<double> > mat = op_.matrix();

  /* fill in with the Laplacian operator plus exp(-x) */
  for (int i=0; i<nLocalRows_; i++)
  {
    int row = lowestLocalRow + i;
    Array<int> colIndices;
    Array<double> colVals;
    if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows_-1))
    {
      colIndices = tuple(row);
      colVals = tuple(1.0);
    }
    else
    {
      colIndices = tuple(row-1, row, row+1);
      colVals = tuple(1.0/h_/h_, 
        -2.0/h_/h_ + exp(-x[i]), 
        1.0/h_/h_);
    }
    mat->addToRow(row, colIndices.size(), 
      &(colIndices[0]), &(colVals[0]));
  }
  Out::os() << tab << "done PBJ::setEvalPoint()" << std::endl;
}
