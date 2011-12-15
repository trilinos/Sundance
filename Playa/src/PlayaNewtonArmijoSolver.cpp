/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaNewtonArmijoSolverImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"


namespace Playa
{

template class NewtonArmijoSolver<double>;

}

#endif
