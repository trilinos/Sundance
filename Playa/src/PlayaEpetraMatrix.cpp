/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  // changed from Impl
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"  // changed from Impl
#include "PlayaIfpackILUOperator.hpp"
#include "PlayaIfpackICCOperator.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaGenericLeftPreconditioner.hpp"
#include "PlayaGenericRightPreconditioner.hpp"

#include "PlayaGenericTwoSidedPreconditioner.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;


EpetraMatrix::EpetraMatrix(const Epetra_CrsGraph& graph,
  const VectorSpace<double>& domain,
  const VectorSpace<double>& range)
  : LinearOpWithSpaces<double>(domain, range),
    matrix_(rcp(new Epetra_CrsMatrix(Copy, graph)))
{}

EpetraMatrix::EpetraMatrix(const RCP<Epetra_CrsMatrix>& mat,
  const VectorSpace<double>& domain,
  const VectorSpace<double>& range)
  : LinearOpWithSpaces<double>(domain, range),
    matrix_(mat)
{}


void EpetraMatrix::apply(
  Teuchos::ETransp applyType,
  const Vector<double>& in,
  Vector<double> out) const
{
  const Epetra_Vector& epIn = EpetraVector::getConcrete(in);
  Epetra_Vector& epOut = EpetraVector::getConcrete(out);
  using Teuchos::rcp_dynamic_cast;

  bool trans = applyType == TRANS;
  
  int ierr = matrix_->Multiply(trans, epIn, epOut);
  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
    "EpetraMatrix::generalApply() detected ierr="
    << ierr << " in matrix_->Multiply()");
}



void EpetraMatrix::addToRow(int globalRowIndex,
  int nElemsToInsert,
  const int* globalColumnIndices,
  const double* elementValues)
{
  int ierr = crsMatrix()->SumIntoGlobalValues(globalRowIndex,
    nElemsToInsert,
    (double*) elementValues,
    (int*) globalColumnIndices);

  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
    "failed to add to row " << globalRowIndex
    << " in EpetraMatrix::addToRow() with nnz="
    << nElemsToInsert 
    << ". Error code was " << ierr);
}



void EpetraMatrix::addToElementBatch(int numRows, 
  int rowBlockSize,
  const int* globalRowIndices,
  int numColumnsPerRow,
  const int* globalColumnIndices,
  const double* values,
  const int* skipRow)
{
  Epetra_CrsMatrix* crs = crsMatrix();

  int numRowBlocks = numRows/rowBlockSize;
  int row = 0;

  for (int rb=0; rb<numRowBlocks; rb++)
  {
    const int* cols = globalColumnIndices + rb*numColumnsPerRow;
    for (int r=0; r<rowBlockSize; r++, row++)
    {
      if (skipRow[row]) continue;
      const double* rowVals = values + row*numColumnsPerRow;
      int ierr=crs->SumIntoGlobalValues(globalRowIndices[row], 
        numColumnsPerRow,
        (double*) rowVals,
        (int*) cols);
      TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, 
        "failed to add to row " << globalRowIndices[row]
        << " in EpetraMatrix::addToRow() with nnz="
        << numColumnsPerRow
        << ". Error code was " << ierr);
    }
  }
}


void EpetraMatrix::zero()
{
  crsMatrix()->PutScalar(0.0);
}



void EpetraMatrix::getILUKPreconditioner(int fillLevels,
  int overlapFill,
  double relaxationValue,
  double relativeThreshold,
  double absoluteThreshold,
  LeftOrRight leftOrRight,
  Preconditioner<double>& rtn) const
{
  RCP<LinearOperatorBase<double> > a = rcp(new IfpackILUOperator(this, 
      fillLevels,
      overlapFill,
      relaxationValue,
      relativeThreshold,
      absoluteThreshold));

  LinearOperator<double> ilu = a;


  if (leftOrRight == Left)
  {
    rtn = new GenericLeftPreconditioner<double>(ilu);
  }
  else
  {
    rtn = new GenericRightPreconditioner<double>(ilu);
  }
}

void EpetraMatrix::getICCPreconditioner(int fillLevels,
  int overlapFill,
  double relaxationValue,
  double relativeThreshold,
  double absoluteThreshold,
  Preconditioner<double>& rtn) const
{
  RCP<LinearOperatorBase<double> > a = rcp(new IfpackICCOperator(this,
      fillLevels,
      overlapFill,
      relaxationValue,
      relativeThreshold,
      absoluteThreshold));

  LinearOperator<double> icc = a;
  LinearOperator<double> iccTrans=icc.transpose();
  //Transpose call and then put in for second argument for below

    rtn = new GenericTwoSidedPreconditioner<double>(icc,iccTrans);
}


void EpetraMatrix::print(std::ostream& os) const 
{
  crsMatrix()->Print(os);
}


string EpetraMatrix::description() const 
{
  std::string rtn = "EpetraMatrix[nRow=" 
    + Teuchos::toString(crsMatrix()->NumGlobalRows())
    + ", nCol=" + Teuchos::toString(crsMatrix()->NumGlobalCols())
    + "]";
  return rtn;
}


Epetra_CrsMatrix* EpetraMatrix::crsMatrix()
{
  return matrix_.get();
}


const Epetra_CrsMatrix* EpetraMatrix::crsMatrix() const 
{
  return matrix_.get();
}


Epetra_CrsMatrix& EpetraMatrix::getConcrete(const LinearOperator<double>& A)
{
  return *Teuchos::dyn_cast<EpetraMatrix>(*A.ptr()).crsMatrix();
}


RCP<const Epetra_CrsMatrix>
EpetraMatrix::getConcretePtr(const LinearOperator<double>& A)
{
  return Teuchos::rcp_dynamic_cast<EpetraMatrix>(A.ptr())->matrix_;
}


void EpetraMatrix::getRow(const int& row, 
  Teuchos::Array<int>& indices, 
  Teuchos::Array<double>& values) const
{
  const Epetra_CrsMatrix* crs = crsMatrix();

  int numEntries;
  int* epIndices;
  double* epValues;

  int info = crs->ExtractGlobalRowView(row, numEntries, epValues, epIndices);
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
    "call to ExtractGlobalRowView not successful");

  indices.resize(numEntries);
  values.resize(numEntries);
  for (int i = 0; i < numEntries; i++)
  {
    indices[i] = *epIndices;
    values[i] = *epValues;
    epIndices++;
    epValues++;
  }
}



