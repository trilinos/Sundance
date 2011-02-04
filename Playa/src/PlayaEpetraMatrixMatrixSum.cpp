/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraMatrixMatrixSum.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaExceptions.hpp"
#include "EpetraExt_MatrixMatrix.h"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif



namespace Playa
{
using namespace Teuchos;


LinearOperator<double> epetraMatrixMatrixSum(
  const LinearOperator<double>& A,
  const LinearOperator<double>& B)
{
  /* Extract the underlying Epetra matrix for A. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);

  /* Extract the underlying Epetra matrix for A. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> B_crs = EpetraMatrix::getConcretePtr(B);
  
  bool transA = false;
  bool transB = false;

  TEST_FOR_EXCEPTION(A.range() != B.range(), RuntimeError,
    "incompatible operand ranges in epetraMatrixMatrixSum()"
    << std::endl << "A.range()=" << A.range()
    << std::endl << "B.range()=" << B.range()
    );
  

  TEST_FOR_EXCEPTION(A.domain() != B.domain(), RuntimeError,
    "incompatible operand domains in epetraMatrixMatrixSum()"
    << std::endl << "A.domain()=" << A.domain()
    << std::endl << "B.domain()=" << B.domain()
    );
  

  /* Get the row map from A. We will need this to build the target matrix C */
  const Epetra_Map* rowmap 
    = transA ? &(A_crs->DomainMap()) : &(A_crs->RowMap());

  /* make the target matrix */
  RCP<Epetra_CrsMatrix> C = rcp(new Epetra_CrsMatrix(Copy, *rowmap, 1));
  Epetra_CrsMatrix* CPtr = C.get();

  /* Carry out the multiplication */
  int ierr 
    = EpetraExt::MatrixMatrix::Add(
      *A_crs, transA, 1.0, 
      *B_crs, transB, 1.0, CPtr);
  TEST_FOR_EXCEPTION(ierr != 0, RuntimeError,
    "EpetraExt Matrix-matrix add failed with error code ierr=" << ierr);

  /* Need to call FillComplete on the result */
  C->FillComplete();

  /* Prepare an operator object for the added matrix */
  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(C, B.domain(), A.range()));
  return rtn;
  
}

}
