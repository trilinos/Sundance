/* @HEADER@ */
//   
 /* @HEADER@ */

#include "PlayaDenseLUSolver.hpp"
#include "PlayaDenseSerialMatrix.hpp"
#include "PlayaLinearOperatorDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaVectorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;

using std::setw;


DenseLUSolver::DenseLUSolver()
  : LinearSolverBase<double>(ParameterList())
{
}

SolverState<double> DenseLUSolver::solve(const LinearOperator<double>& op,
  const Vector<double>& rhs,
  Vector<double>& soln) const
{
  return denseSolve(op, rhs, soln);
}
