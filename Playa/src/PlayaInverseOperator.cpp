/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaInverseOperatorImpl.hpp"

namespace Playa
{
template class InverseOperator<double>;

template
LinearOperator<double> 
inverse(const LinearOperator<double>& op, 
  const LinearSolver<double>& solver);
}

#endif
