/* @HEADER@ */
//   
 /* @HEADER@ */


#ifndef RANDOMSPARSEMATRIX_BUILDER_IMPL_HPP
#define RANDOMSPARSEMATRIX_BUILDER_IMPL_HPP

#include "PlayaRandomSparseMatrixBuilderDecl.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"
#include "PlayaLoadableMatrix.hpp"


using namespace Playa;
using namespace Teuchos;


namespace Playa
{

template <class Scalar> 
inline RandomSparseMatrixBuilder<Scalar>
::RandomSparseMatrixBuilder(int nLocalRows, int nLocalCols,
  double onProcDensity,
  double offProcDensity,
  const VectorType<double>& type)
  : OperatorBuilder<double>(nLocalRows, nLocalCols, type), op_()
{
  initOp(onProcDensity, offProcDensity);
}


template <class Scalar> 
inline RandomSparseMatrixBuilder<Scalar>
::RandomSparseMatrixBuilder(const VectorSpace<Scalar>& d,
  const VectorSpace<Scalar>& r,
  double onProcDensity,
  double offProcDensity,
  const VectorType<double>& type)
  : OperatorBuilder<double>(d, r, type), op_()
{
  initOp(onProcDensity, offProcDensity);
}


template <class Scalar> 
inline void RandomSparseMatrixBuilder<Scalar>
::initOp(double onProcDensity,
  double offProcDensity)
{
  int rank = MPIComm::world().getRank();
  int nProc = MPIComm::world().getNProc();

  RCP<MatrixFactory<double> > mFact 
    = this->vecType().createMatrixFactory(this->domain(), this->range());

  int colDimension = this->domain().dim();
  int rowDimension = this->range().dim();
  int numLocalCols = colDimension / nProc;
  int numLocalRows = rowDimension / nProc;
  int lowestLocalRow = numLocalRows * rank;

  int lowestLocalCol = numLocalCols * rank;
  int highestLocalCol = numLocalCols * (rank+1) - 1;


  IncrementallyConfigurableMatrixFactory* icmf 
    = dynamic_cast<IncrementallyConfigurableMatrixFactory*>(mFact.get());
  Array<Array<int> > colIndices(numLocalRows);
  for (int i=0; i<numLocalRows; i++)
  {
    int row = lowestLocalRow + i;

    Array<int>& cols = colIndices[i];

    while (cols.size() == 0)
    {
      for (int j=0; j<colDimension; j++)
      {
        double acceptProb;
        if (j >= lowestLocalCol && j <= highestLocalCol)
        {
          acceptProb = onProcDensity;
        }
        else
        {
          acceptProb = offProcDensity;
        }
        double p = 0.5*(ScalarTraits<double>::random() + 1.0);

        if (p < acceptProb)
        {
          cols.append(j);
        }
      }
      if (cols.size()>0)
      {
        icmf->initializeNonzerosInRow(row, colIndices[i].size(),
          &(colIndices[i][0]));
      }
    }
        
  }
  icmf->finalize();
      
  op_ = mFact->createMatrix();
      
  RCP<LoadableMatrix<double> > mat = op_.matrix();

  /* fill in with the Laplacian operator */
  for (int i=0; i<numLocalRows; i++)
  {
    int row = lowestLocalRow + i;
    const Array<int>& cols = colIndices[i];
    Array<Scalar> colVals(cols.size());
    for (int j=0; j<cols.size(); j++)
    {
      colVals[j] = ScalarTraits<Scalar>::random();
    }
    if (cols.size() > 0)
    {
      mat->addToRow(row, colIndices[i].size(), 
        &(colIndices[i][0]), &(colVals[0]));
    }
  }
}
}

#endif
