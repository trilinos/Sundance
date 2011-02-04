/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaDefaultBlockVectorSpaceImpl.hpp"

template class Playa::DefaultBlockVectorSpace<double>;

namespace Playa
{


template VectorSpace<double> blockSpace(
  const VectorSpace<double>& v1,
  const VectorSpace<double>& v2);


template VectorSpace<double> blockSpace(
  const VectorSpace<double>& v1,
  const VectorSpace<double>& v2,
  const VectorSpace<double>& v3);

template VectorSpace<double> blockSpace(
  const VectorSpace<double>& v1,
  const VectorSpace<double>& v2,
  const VectorSpace<double>& v3,
  const VectorSpace<double>& v4);

template VectorSpace<double> blockSpace(const Array<VectorSpace<double> >& x);

}


#endif
