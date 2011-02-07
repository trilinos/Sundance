/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTOROPSDECL_HPP
#define PLAYA_VECTOROPSDECL_HPP

#include "PlayaDefs.hpp"

namespace Playa
{

template <class Scalar> class Vector;

/** */
template <class Scalar> 
Scalar minloc(const Vector<Scalar>& x, int& gni);

/** */
template <class Scalar> 
Scalar maxloc(const Vector<Scalar>& x, int& gni);

/** */
template <class Scalar> 
Scalar minlocWithBound(const Scalar& lowerBound, 
  const Vector<Scalar>& x, int& gni);

/** */
template <class Scalar> 
Scalar maxlocWithBound(const Scalar& upperBound, 
  const Vector<Scalar>& x, int& gni);

}

 

#endif
