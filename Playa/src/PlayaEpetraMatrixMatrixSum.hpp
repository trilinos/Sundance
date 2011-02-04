/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_EPETRA_MATRIXMATRIXSUM_HPP
#define PLAYA_EPETRA_MATRIXMATRIXSUM_HPP


#include "PlayaLinearOperatorDecl.hpp"

namespace Playa
{

/** \relates LinearOperator */
LinearOperator<double> epetraMatrixMatrixSum(
  const LinearOperator<double>& A,
  const LinearOperator<double>& B);


}

#endif
