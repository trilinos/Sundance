#include "PlayaOptState.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif


namespace Playa
{


OptState::OptState(const Vector<double>& xCur,
  const double& fCur,
  const Vector<double>& gradCur)
  : status_(Opt_Continue),
    iter_(0),
    xCur_(xCur),
    xPrev_(),
    gradCur_(gradCur),
    gradPrev_(),
    fCur_(fCur),
    fPrev_(HUGE_VAL)
{}

void OptState::update(const Vector<double>& xNew, 
  const Vector<double>& gradNew, 
  const double& fNew)
{
  xPrev_.acceptCopyOf(xCur_);
  gradPrev_.acceptCopyOf(gradCur_);
  fPrev_ = fCur_;
  
  xCur_.acceptCopyOf(xNew);
  gradCur_.acceptCopyOf(gradNew);
  fCur_ = fNew;
  
  iter_++;
}



}
