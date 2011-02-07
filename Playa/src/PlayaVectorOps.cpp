/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaVectorOpsImpl.hpp"

namespace Playa
{

/** */
template double minloc(const Vector<double>& x, int& gni);

/** */
template double maxloc(const Vector<double>& x, int& gni);

/** */
template double minlocWithBound(const double& lowerBound, 
  const Vector<double>& x, int& gni);

/** */
template double maxlocWithBound(const double& upperBound, 
  const Vector<double>& x, int& gni);
  

}

#endif

