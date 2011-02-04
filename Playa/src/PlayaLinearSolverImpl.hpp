/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_LINEARSOLVERIMPL_HPP
#define PLAYA_LINEARSOLVERIMPL_HPP

#include "PlayaLinearSolverDecl.hpp"
#include "PlayaLinearSolverBaseImpl.hpp"
#include "PlayaPreconditionerFactory.hpp"

using namespace Playa;
using namespace Teuchos;



template <class Scalar> inline
void LinearSolver<Scalar>::setUserPrec(const PreconditionerFactory<Scalar>& pf)
{
  this->ptr()->setUserPrec(pf);
}

template <class Scalar> inline
void LinearSolver<Scalar>::setUserPrec(const LinearOperator<Scalar>& P,
  const LinearSolver<Scalar>& pSolver)
{
  this->ptr()->setUserPrec(P, pSolver);
}

#endif
