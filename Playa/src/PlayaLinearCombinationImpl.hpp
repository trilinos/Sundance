/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LINEARCOMBINATIONIMPL_HPP
#define PLAYA_LINEARCOMBINATIONIMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif


namespace Playa
{

using Playa::Out;
using Playa::Tabs;
using std::endl;

/*=========================================================================
 * Copy constructors
 * ======================================================================== */


template <class Scalar> 
template <int N>
inline
Vector<Scalar>::Vector(const LCN<Scalar, N>& x)
  : Playa::Handle<VectorBase<Scalar> >()
{
  this->ptr() = x.eval().ptr();
}

/*=========================================================================
 * Assignment operators
 * ======================================================================== */

template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator=(const LCN<Scalar, 1>& lc)
{
  const Vector<Scalar>& other = lc.vec(0);
  const Scalar& alpha = lc.coeff(0);

  if (this->ptr().get() == other.ptr().get())
  {
    /* If LHS and RHS vectors are identical, simply multiply LHS
     * by the scalar constant */
    if (alpha == 0.0) this->zero();
    else if (alpha != 1.0) this->scale(alpha);
  }
  else if (this->ptr().get() != 0 && this->space() == other.space())
  {
    /* If the vectors are distinct but from the same space, use the
     * update operation to compute (*this) = 0.0*(*this) + alpha*other */ 
    this->ptr()->update(alpha, other.ptr().get(), 0.0);
  }
  else
  {
    /* If the vectors are from different spaces, or if the LHS is null,
     * copy the RHS vector into the LHS and scale */
    Vector<Scalar> cp = other.copy();
    this->ptr() = cp.ptr();
    this->scale(alpha);
    
  }
  return *this;
}  


template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator=(const LCN<Scalar, 2>& lc)
{
  const Vector<Scalar>& x = lc.vec(0);
  const Scalar& alpha = lc.coeff(0);
  const Vector<Scalar>& y = lc.vec(1);
  const Scalar& beta = lc.coeff(1);

  TEST_FOR_EXCEPTION(!y.space().isCompatible(x.space()),
    std::runtime_error,
    "Spaces x=" << x.space() << " and y="
    << y.space() << " are not compatible in operator=(a*x + b*y)");

  if (this->ptr().get() != 0 &&  this->space() == x.space())
  {
    /* If the LHS exists and is from the same space as the RHS vectors, use the
     * update operation to compute 
     *      (*this) = 0.0*(*this) + alpha*x + beta*y 
     */ 
    this->update(alpha, x, beta, y, 0.0);
  }
  else
  {
    /* If the vectors are from different spaces, or if the LHS is null,
     * form the RHS vector and overwrite the LHS's ptr with it. */
    Vector e = lc.eval();
    this->ptr() = e.ptr();
  }

  return *this;
}  

template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator=(const LCN<Scalar, 3>& lc)
{
  const Vector<Scalar>& x = lc.vec(0);
  const Scalar& alpha = lc.coeff(0);
  const Vector<Scalar>& y = lc.vec(1);
  const Scalar& beta = lc.coeff(1);
  const Vector<Scalar>& z = lc.vec(2);
  const Scalar& gamma = lc.coeff(2);

  TEST_FOR_EXCEPTION(!y.space().isCompatible(x.space()),
    std::runtime_error,
    "Spaces x=" << x.space() << " and y="
    << y.space() << " are not compatible in operator=(a*x + b*y + c*z)");

  TEST_FOR_EXCEPTION(!z.space().isCompatible(x.space()),
    std::runtime_error,
    "Spaces x=" << x.space() << " and z="
    << z.space() << " are not compatible in operator=(a*x + b*y + c*z)");

  if (this->ptr().get() != 0 &&  this->space() == x.space())
  {
    /* If the LHS exists and is from the same space as the RHS vectors, use the
     * update operation to compute 
     *      (*this) = 0.0*(*this) + alpha*x + beta*y + c*z
     */ 
    this->update(alpha, x, beta, y, gamma, z, 0.0);
  }
  else
  {
    /* If the vectors are from different spaces, or if the LHS is null,
     * form the RHS vector and overwrite the LHS's ptr with it. */
    Vector e = lc.eval();
    this->ptr() = e.ptr();
  }

  return *this;
}  

template <class Scalar> 
template <int N> inline
Vector<Scalar>& Vector<Scalar>::operator=(const LCN<Scalar, N>& lc)
{
  Vector e = lc.eval();
  this->ptr() = e.ptr();
    
  return *this;
}  

  
/*=========================================================================
 * Reflexive addition operators
 * ======================================================================== */



template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator+=(const LCN<Scalar, 1>& lc)
{
  const Vector<Scalar>& other = lc.vec(0);
  const Scalar& alpha = lc.coeff(0);

  TEST_FOR_EXCEPTION(!this->space().isCompatible(other.space()),
    std::runtime_error,
    "Spaces this=" << this->space() << " and other="
    << other.space() << " are not compatible in operator+=()");

  this->ptr()->update(alpha, other.ptr().get(), 1.0);

  return *this;
}  


template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator+=(const LCN<Scalar, 2>& lc)
{
  const Vector<Scalar>& x = lc.vec(0);
  const Scalar& alpha = lc.coeff(0);
  const Vector<Scalar>& y = lc.vec(1);
  const Scalar& beta = lc.coeff(1);

  TEST_FOR_EXCEPTION(!this->space().isCompatible(x.space()),
    std::runtime_error,
    "Spaces this=" << this->space() << " and other="
    << x.space() << " are not compatible in operator+=()");

  this->update(alpha, x, beta, y, 1.0);

  return *this;
}  


template <class Scalar> 
template <int N> inline 
Vector<Scalar>& Vector<Scalar>::operator+=(const LCN<Scalar, N>& lc)
{
  Vector<Scalar> e = lc.eval();
  *this += e;
  return *this;
}  

  

/*=========================================================================
 * Reflexive subtraction operators
 * ======================================================================== */



template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator-=(const LCN<Scalar, 1>& lc)
{
  const Vector<Scalar>& other = lc.vec(0);
  const Scalar& alpha = -lc.coeff(0);

  TEST_FOR_EXCEPTION(!this->space().isCompatible(other.space()),
    std::runtime_error,
    "Spaces this=" << this->space() << " and other="
    << other.space() << " are not compatible in operator+=()");

  this->ptr()->update(alpha, other.ptr().get(), 1.0);

  return *this;
}  


template <class Scalar> inline
Vector<Scalar>& Vector<Scalar>::operator-=(const LCN<Scalar, 2>& lc)
{

  const Vector<Scalar>& x = lc.vec(0);
  const Scalar& alpha = -lc.coeff(0);
  const Vector<Scalar>& y = lc.vec(1);
  const Scalar& beta = -lc.coeff(1);

  TEST_FOR_EXCEPTION(!this->space().isCompatible(x.space()),
    std::runtime_error,
    "Spaces this=" << this->space() << " and other="
    << x.space() << " are not compatible in operator+=()");

  this->update(alpha, x, beta, y, 1.0);

  return *this;
}  


template <class Scalar> 
template <int N> inline 
Vector<Scalar>& Vector<Scalar>::operator-=(const LCN<Scalar, N>& lc)
{
  Vector<Scalar> e = lc.eval();
  *this -= e;

  return *this;
}  

/*======================================================================
 *
 *    operator times vector or LC
 *
 *======================================================================*/

/* */
template <class Scalar> 
Vector<Scalar> operator*(const LinearOperator<Scalar>& A,
  const Vector<Scalar>& x)
{
  Vector<Scalar> rtn;
  A.apply(x, rtn);
  return rtn;
}


/* */
template <class Scalar, int N> 
Vector<Scalar> operator*(const LinearOperator<Scalar>& A,
  const LCN<Scalar, N>& x)
{
  return A*x.eval();
}
  
/* */
template <class Scalar> 
Vector<Scalar> operator*(const LinearOperator<Scalar>& A,
  const LCN<Scalar, 1>& x)
{
  Vector<Scalar> rtn = A*x.vec(0);
  if (x.coeff(0)!=1.0) rtn.scale(x.coeff(0));
  return rtn;
}


/*======================================================================
 *
 *    scalar times vector
 *
 *======================================================================*/

/* \relates LCN scalar * vec */
template <class Scalar> inline
LCN<Scalar, 1> operator*(const Scalar& alpha, 
  const Vector<Scalar>& x)
{
  return LCN<Scalar, 1>(alpha, x);
}

/* \relates LCN vec * scalar */
template <class Scalar> inline
LCN<Scalar, 1> operator*(const Vector<Scalar>& x, 
  const Scalar& alpha)
{
  return LCN<Scalar, 1>(alpha, x);
}


/* \relates LCN vec / scalar */
template <class Scalar> inline
LCN<Scalar, 1> operator/(const Vector<Scalar>& x, 
  const Scalar& alpha)
{
  return LCN<Scalar, 1>(1.0/alpha, x);
}


/*======================================================================
 *
 *    scalar times LC
 *
 *======================================================================*/


/* */
template <class Scalar, int N> inline
LCN<Scalar, N> operator*(const LCN<Scalar, N>& lc, const Scalar& beta)
{
  LCN<Scalar, N> rtn(lc);
  rtn.multiply(beta);
  return rtn;
}

/* */
template <class Scalar, int N> inline
LCN<Scalar, N> operator*(const Scalar& beta, const LCN<Scalar, N>& lc)
{
  LCN<Scalar, N> rtn(lc);
  rtn.multiply(beta);
  return rtn;
}

/* */
template <class Scalar, int N> inline
LCN<Scalar, N> operator/(const LCN<Scalar, N>& lc, const Scalar& beta)
{
  LCN<Scalar, N> rtn(lc);
  rtn.multiply(1.0/beta);
  return rtn;
}


  
/*======================================================================
 *
 *  vector plus/minus vector
 *
 *======================================================================*/


/* */
template <class Scalar> inline
LCN<Scalar, 2> operator+(const Vector<Scalar>& x, 
  const Vector<Scalar>& y)
{
  return LCN<Scalar, 2>(1.0, x, 1.0, y);
}

/* */
template <class Scalar> inline
LCN<Scalar, 2> operator-(const Vector<Scalar>& x, 
  const Vector<Scalar>& y)
{
  return LCN<Scalar, 2>(1.0, x, -1.0, y);
}


/*======================================================================
 *
 *  LC plus LC
 *
 *======================================================================*/


/* */
template <class Scalar, int N, int M> inline
LCN<Scalar, N+M> operator+(const LCN<Scalar, N>& f, const LCN<Scalar, M>& g)
{
  LCN<Scalar, N+M> rtn;
  for (int i=0; i<N; i++) rtn.set(i, f.coeff(i), f.vec(i));
  for (int i=0; i<M; i++) rtn.set(i+N, g.coeff(i), g.vec(i));
  return rtn;
}

/* */
template <class Scalar, int N> inline
LCN<Scalar, N+1> operator+(const LCN<Scalar, N>& f, const Vector<Scalar>& g)
{
  LCN<Scalar, N+1> rtn;
  for (int i=0; i<N; i++) rtn.set(i, f.coeff(i), f.vec(i));
  rtn.set(N, 1.0, g);
  return rtn;
}

/* */
template <class Scalar, int N> inline
LCN<Scalar, N+1> operator+(const Vector<Scalar>& f, const LCN<Scalar, N>& g)
{
  LCN<Scalar, N+1> rtn;
  rtn.set(0, 1.0, f);
  for (int i=0; i<N; i++) rtn.set(i+1, g.coeff(i), g.vec(i));

  return rtn;
}

/* */

template <class Scalar> inline
LCN<Scalar, 2> operator+(const LCN<Scalar, 1>& lc, const Vector<Scalar>& x)
{
  return LCN<Scalar, 2>(lc, 1.0, x);
}


/* */

template <class Scalar> inline
LCN<Scalar, 2> operator+(const Vector<Scalar>& x, const LCN<Scalar, 1>& lc)
{
  return LCN<Scalar, 2>(1.0, x, lc);
}

/* */

template <class Scalar> inline
LCN<Scalar, 2> operator+(const LCN<Scalar, 1>& ax, const LCN<Scalar, 1>& by)
{
  return LCN<Scalar, 2>(ax, by);
}


/* */

template <class Scalar> inline
LCN<Scalar, 3> operator+(const LCN<Scalar, 1>& ax, const LCN<Scalar, 2>& bycz)
{
  return LCN<Scalar, 3>(ax, bycz);
}

/* */

template <class Scalar> inline
LCN<Scalar, 3> operator+(const Vector<Scalar>& x, const LCN<Scalar, 2>& bycz)
{
  return LCN<Scalar, 3>(LCN<Scalar, 1>(1.0, x), bycz);
}


/* */

template <class Scalar> inline
LCN<Scalar, 3> operator+(const LCN<Scalar, 2>& axby, const LCN<Scalar, 1>& cz)
{
  return LCN<Scalar, 3>(axby, cz);
}


/* */

template <class Scalar> inline
LCN<Scalar, 3> operator+(const LCN<Scalar, 2>& axby, const Vector<Scalar>& z)
{
  return LCN<Scalar, 3>(axby, LCN<Scalar,1>(1.0,z));
}



/*======================================================================
 *
 *  LC minus LC
 *
 *======================================================================*/


/* */
template <class Scalar, int N, int M> inline
LCN<Scalar, N+M> operator-(const LCN<Scalar, N>& f, const LCN<Scalar, M>& g)
{
  LCN<Scalar, N+M> rtn;
  for (int i=0; i<N; i++) rtn.set(i, f.coeff(i), f.vec(i));
  for (int i=0; i<M; i++) rtn.set(i+N, -g.coeff(i), g.vec(i));
  return rtn;
}

/* */
template <class Scalar, int N> inline
LCN<Scalar, N+1> operator-(const LCN<Scalar, N>& f, const Vector<Scalar>& g)
{
  LCN<Scalar, N+1> rtn;
  for (int i=0; i<N; i++) rtn.set(i, f.coeff(i), f.vec(i));
  rtn.set(N, -1.0, g);
  return rtn;
}

/* */
template <class Scalar, int N> inline
LCN<Scalar, N+1> operator-(const Vector<Scalar>& f, const LCN<Scalar, N>& g)
{
  LCN<Scalar, N+1> rtn;
  rtn.set(0, 1.0, f);
  for (int i=0; i<N; i++) rtn.set(i+1, -g.coeff(i), g.vec(i));

  return rtn;
}

/* */

template <class Scalar> inline
LCN<Scalar, 2> operator-(const LCN<Scalar, 1>& lc, const Vector<Scalar>& x)
{
  return LCN<Scalar, 2>(lc, -1.0, x);
}


/* */

template <class Scalar> inline
LCN<Scalar, 2> operator-(const Vector<Scalar>& x, const LCN<Scalar, 1>& lc)
{
  return LCN<Scalar, 2>(1.0, x, -lc.coeff(0), lc.vec(0));
}

/* */

template <class Scalar> inline
LCN<Scalar, 2> operator-(const LCN<Scalar, 1>& ax, const LCN<Scalar, 1>& by)
{
  return LCN<Scalar, 2>(ax, -by);
}


/* */

template <class Scalar> inline
LCN<Scalar, 3> operator-(const LCN<Scalar, 1>& ax, const LCN<Scalar, 2>& bycz)
{
  return LCN<Scalar, 3>(ax, -bycz);
}

/* */

template <class Scalar> inline
LCN<Scalar, 3> operator-(const Vector<Scalar>& x, const LCN<Scalar, 2>& bycz)
{
  return LCN<Scalar, 3>(LCN<Scalar, 1>(1.0, x), -bycz);
}


/* */

template <class Scalar> inline
LCN<Scalar, 3> operator-(const LCN<Scalar, 2>& axby, const LCN<Scalar, 1>& cz)
{
  return LCN<Scalar, 3>(axby, -cz);
}


/* */

template <class Scalar> inline
LCN<Scalar, 3> operator-(const LCN<Scalar, 2>& axby, const Vector<Scalar>& z)
{
  return LCN<Scalar, 3>(axby, LCN<Scalar,1>(-1.0,z));
}




/*======================================================================
 *
 *  Operations on LC
 *
 *======================================================================*/

/* */
template <class Scalar, int N> inline
Scalar norm1(const LCN<Scalar, N>& lc) 
{
  return lc.norm1();
}

/* */
template <class Scalar, int N> inline
Scalar norm2(const LCN<Scalar, N>& lc) 
{
  return lc.norm2();
}

/* */
template <class Scalar, int N> inline
Scalar normInf(const LCN<Scalar, N>& lc) 
{
  return lc.normInf();
}

/* */
template <class Scalar, int N> inline
Vector<Scalar> abs(const LCN<Scalar, N>& lc) 
{
  return lc.abs();
}

/* */
template <class Scalar, int N> inline
Scalar min(const LCN<Scalar, N>& lc) 
{
  return lc.min();
}

/* */
template <class Scalar, int N> inline
Scalar max(const LCN<Scalar, N>& lc) 
{
  return lc.max();
}

/* */
template <class Scalar, int N> inline
Vector<Scalar> reciprocal(const LCN<Scalar, N>& lc) 
{
  return lc.reciprocal();
}


  
}



#endif
