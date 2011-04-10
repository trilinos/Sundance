/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTOROPSDECL_HPP
#define PLAYA_VECTOROPSDECL_HPP

#include "PlayaDefs.hpp"

namespace Playa
{

template <class Scalar> class Vector;

/** \relates Vector \brief Return minimum element and its location */
template <class Scalar> 
Scalar minloc(const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Return maximum element and its location*/
template <class Scalar> 
Scalar maxloc(const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Return minimum element greater than a specified
bound, and its location*/
template <class Scalar> 
Scalar minlocWithBound(const Scalar& lowerBound, 
  const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Return maximum element less than a specified
bound, and its location*/
template <class Scalar> 
Scalar maxlocWithBound(const Scalar& upperBound, 
  const Vector<Scalar>& x, int& gni);

/** \relates Vector \brief Compute the Euclidean norm of a vector */
template <class Scalar>
Scalar norm2(const Vector<Scalar>& x);

/** \relates Vector \brief Compute the one-norm of a vector */
template <class Scalar>
Scalar norm1(const Vector<Scalar>& x);

/** \relates Vector \brief Compute the infinity norm of a vector */
template <class Scalar>
Scalar normInf(const Vector<Scalar>& x);

/** \relates Vector 
 * \brief Compute the Euclidean distance between two vectors */
template <class Scalar>
Scalar norm2Dist(const Vector<Scalar>& x, const Vector<Scalar>& y);

/** \relates Vector \brief Compute the one-norm distance between two vectors */
template <class Scalar>
Scalar norm1Dist(const Vector<Scalar>& x, const Vector<Scalar>& y);

/** \relates Vector
 *  \brief Compute the infinity-norm distance between two vectors */
template <class Scalar>
Scalar normInfDist(const Vector<Scalar>& x, const Vector<Scalar>& y);

}

 

#endif
