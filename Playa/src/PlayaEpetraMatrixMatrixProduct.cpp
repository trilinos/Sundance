/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraMatrixMatrixProduct.hpp"
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


LinearOperator<double> epetraLeftScale(
  const Vector<double>& d,
  const LinearOperator<double>& A)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);
  
  /* Make a deep copy of A */
  RCP<Epetra_CrsMatrix> mtxCopy = rcp(new Epetra_CrsMatrix(*A_crs));

  /* Extract the underlying Epetra vector. Type checking is done
   * internally, so we need no error check here. */
  const Epetra_Vector& epv = EpetraVector::getConcrete(d);
  
  /* Scale the copy */
  mtxCopy->LeftScale(epv);

  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(mtxCopy, A.domain(), A.range()));
  return rtn;
  
}

LinearOperator<double> epetraRightScale(
  const LinearOperator<double>& A,
  const Vector<double>& d)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);
  
  /* Make a deep copy of A */
  RCP<Epetra_CrsMatrix> mtxCopy = rcp(new Epetra_CrsMatrix(*A_crs));

  /* Extract the underlying Epetra vector. Type checking is done
   * internally, so we need no error check here. */
  const Epetra_Vector& epv = EpetraVector::getConcrete(d);
  
  /* Scale the copy */
  mtxCopy->RightScale(epv);

  /* Prepare an operator object for the scaled matrix */
  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(mtxCopy, A.domain(), A.range()));
  return rtn;
  
}


LinearOperator<double> epetraMatrixMatrixProduct(
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
  

  /* Get the row map from A. We will need this to build the target matrix C */
  const Epetra_Map* rowmap 
    = transA ? &(A_crs->DomainMap()) : &(A_crs->RowMap());

  /* make the target matrix */
  RCP<Epetra_CrsMatrix> C = rcp(new Epetra_CrsMatrix(Copy, *rowmap, 1));

  /* Carry out the multiplication */
  int ierr 
    = EpetraExt::MatrixMatrix::Multiply(*A_crs, transA, *B_crs, transB, *C);
  TEST_FOR_EXCEPTION(ierr != 0, RuntimeError,
    "EpetraExt Matrix-matrix multiply failed with error code ierr=" << ierr);

  /* Prepare an operator object for the scaled matrix */
  RCP<LinearOperatorBase<double> > rtn 
    = rcp(new EpetraMatrix(C, B.domain(), A.range()));
  return rtn;
  
}

}
