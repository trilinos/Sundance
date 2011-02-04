/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaLinearSolverBaseImpl.hpp"


template class Playa::LinearSolverBase<double>;
template void Playa::LinearSolverBase<double>::setParameter(const Teuchos::ParameterList& pl, int* val, const std::string& name);
template void Playa::LinearSolverBase<double>::setParameter(const Teuchos::ParameterList& pl, bool* val, const std::string& name);
template void Playa::LinearSolverBase<double>::setParameter(const Teuchos::ParameterList& pl, double* val, const std::string& name);

#endif
