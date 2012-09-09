/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_LINEARSOLVERBASEIMPL_HPP
#define PLAYA_LINEARSOLVERBASEIMPL_HPP

#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaPreconditioner.hpp"
#include "PlayaPreconditionerFactory.hpp"
#include "Teuchos_ParameterList.hpp"


using namespace Teuchos;


namespace Playa
{

template <class Scalar> inline
const ParameterList& LinearSolverBase<Scalar>::parameters() const 
{return params_;}

template <class Scalar> inline
LinearSolverBase<Scalar>::LinearSolverBase(const ParameterList& params)
  : ObjectWithVerbosity(), params_(params) 
{
  if (this->parameters().isParameter(this->verbosityParam()))
  {
    this->setVerb(this->parameters().template get<int>(this->verbosityParam()));
  }
}

template <class Scalar> inline
ParameterList& LinearSolverBase<Scalar>::parameters() {return params_;}


template <class Scalar> inline
string LinearSolverBase<Scalar>::verbosityParam() const {return "Verbosity";}

template <class Scalar>
template <typename T> inline
void LinearSolverBase<Scalar>::setParameter(const ParameterList& params,
  T* dataPtr,
  const std::string& name)
{
  if (!params.isParameter(name)) return;

  TEUCHOS_TEST_FOR_EXCEPTION(!params.template isType<T>(name), std::runtime_error,
    "invalid type for parameter [" << name << "]"); 

  *dataPtr = params.template get<T>(name);
}

template <class Scalar> inline
void LinearSolverBase<Scalar>::setUserPrec(const PreconditionerFactory<Scalar>& pf)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "User-defined preconditioning not allowed for generic "
    "linear solver subtypes");
}

template <class Scalar> inline
void LinearSolverBase<Scalar>::setUserPrec(const LinearOperator<Scalar>& P,
  const LinearSolver<Scalar>& pSolver)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
    "User-defined preconditioning not allowed for generic "
    "linear solver subtypes");
}

}

#endif
