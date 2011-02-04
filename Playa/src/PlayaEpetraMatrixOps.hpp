/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_EPETRAMATRIXOPS_HPP
#define PLAYA_EPETRAMATRIXOPS_HPP

#include "PlayaLinearOperatorDecl.hpp"

namespace Playa
{

/** \relates EpetraMatrix */
Vector<double> getEpetraDiagonal(const LinearOperator<double>& A);

/** \relates EpetraMatrix */
LinearOperator<double> makeEpetraDiagonalMatrix(const Vector<double>& d);
  

}

#endif
