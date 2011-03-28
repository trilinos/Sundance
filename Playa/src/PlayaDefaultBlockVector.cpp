/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaDefaultBlockVectorImpl.hpp"

template class Playa::DefaultBlockVector<double>;

namespace Playa
{
template Vector<double> blockVector(
  const Vector<double>& v1);

template Vector<double> blockVector(
  const Vector<double>& v1,
  const Vector<double>& v2);


template Vector<double> blockVector(
  const Vector<double>& v1,
  const Vector<double>& v2,
  const Vector<double>& v3);

template Vector<double> blockVector(
  const Vector<double>& v1,
  const Vector<double>& v2,
  const Vector<double>& v3,
  const Vector<double>& v4);

template Vector<double> blockVector(const Array<Vector<double> >& x);
}

#endif
