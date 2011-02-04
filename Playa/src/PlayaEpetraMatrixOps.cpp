/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraMatrixFactory.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"  // changed from Impl
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"  // changed from Impl

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;




namespace Playa
{

Vector<double> getEpetraDiagonal(const LinearOperator<double>& A)
{
  /* Extract the underlying Epetra matrix. Type checking is done
   * ny rcp_dynamic_cast, so we need no error check here. */
  RCP<const Epetra_CrsMatrix> A_crs = EpetraMatrix::getConcretePtr(A);

  VectorSpace<double> rowSpace = A.domain();
  Vector<double> rtn = rowSpace.createMember();

  Epetra_Vector* xPtr = EpetraVector::getConcretePtr(rtn);
  A_crs->ExtractDiagonalCopy(*xPtr);

  return rtn;
}


LinearOperator<double> makeEpetraDiagonalMatrix(const Vector<double>& d)
{
  VectorSpace<double> space = d.space();
  RCP<const EpetraVectorSpace> eps 
    = rcp_dynamic_cast<const EpetraVectorSpace>(space.ptr());

  EpetraMatrixFactory mf(eps, eps);

  int nLocal = space.numLocalElements();
  int offset = space.baseGlobalNaturalIndex();
  for (int i=0; i<nLocal; i++)
  {
    int rowIndex = offset + i;
    mf.initializeNonzerosInRow(rowIndex, 1, &rowIndex);
  }

  mf.finalize();
  LinearOperator<double> rtn = mf.createMatrix();

  RCP<EpetraMatrix> epm = rcp_dynamic_cast<EpetraMatrix>(rtn.ptr());
  epm->zero();
  
  for (int i=0; i<nLocal; i++)
  {
    int rowIndex = offset + i;
    double val = d[i];
    epm->addToRow(rowIndex, 1, &rowIndex, &val);
  }
  
  return rtn;
}

}
