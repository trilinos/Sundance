/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaHeatOperator1D.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;


HeatOperator1D::HeatOperator1D(int nLocalRows, const VectorType<double>& type)
  : OperatorBuilder<double>(nLocalRows, type), 
    matrixFactory_(),
    cj_(0.0),
    h_(0.0),
    nLocalRows_(nLocalRows),
    lowestLocalRow_(0)
{
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  matrixFactory_ = vecType().createMatrixFactory(domain(), range());

  int n = domain().dim();
  h_ = 1.0/(n - 1.0);

  lowestLocalRow_ = nLocalRows * rank;

  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(matrixFactory_.get());
  for (int i=0; i<nLocalRows; i++)
    {
      int row = lowestLocalRow_ + i;
      Array<int> colIndices;
      if ((rank==0 && i==0) || (rank==nProc-1 && i==nLocalRows-1))
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
}

 

LinearOperator<double> HeatOperator1D::getOp() const  
{

  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();
  LinearOperator<double> rtn = matrixFactory_->createMatrix();
      
  RCP<LoadableMatrix<double> > mat = rtn.matrix();

  /* fill in with the heat operator given the current value of c_j */
  for (int i=0; i<nLocalRows_; i++)
    {
      int row = lowestLocalRow_ + i;
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
          colVals = tuple(-1.0/h_/h_, 2.0/h_/h_ + cj_, -1.0/h_/h_);
        }
      mat->addToRow(row, colIndices.size(), 
                    &(colIndices[0]), &(colVals[0]));
    }

  return rtn;
}
