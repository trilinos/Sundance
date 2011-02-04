/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_EPETRA_MATRIXMATRIXPRODUCT_HPP
#define PLAYA_EPETRA_MATRIXMATRIXPRODUCT_HPP


#include "PlayaLinearOperatorDecl.hpp"

namespace Playa
{

/** \relates LinearOperator */
LinearOperator<double> epetraMatrixMatrixProduct(
  const LinearOperator<double>& A,
  const LinearOperator<double>& B);

/** \relates LinearOperator */
LinearOperator<double> epetraLeftScale(
  const Vector<double>& d,
  const LinearOperator<double>& A);

/** \relates LinearOperator */
LinearOperator<double> epetraRightScale(
  const LinearOperator<double>& A,
  const Vector<double>& d);

}

#endif
