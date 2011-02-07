/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaVectorImpl.hpp"

namespace Playa
{

template class Vector<double>;

template LoadableVector<double>* loadable(Vector<double> vec);

template 
double* dataPtr(Vector<double> vec) ;

template 
const double* dataPtr(const Vector<double>& vec) ;
}

#endif
