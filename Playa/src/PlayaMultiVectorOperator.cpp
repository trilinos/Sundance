/* @HEADER@ */

/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaMultiVectorOperatorImpl.hpp"

namespace Playa
{
template class MultiVectorOperator<double>;

template 
LinearOperator<double> multiVectorOperator(
  const Teuchos::Array<Vector<double> >& cols,
  const VectorSpace<double>& domain);

}

#endif

