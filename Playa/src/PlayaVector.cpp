/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaVectorImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"

namespace Playa
{

template class Vector<double>;

template LoadableVector<double>* loadable(Vector<double> vec);

template 
double* dataPtr(Vector<double> vec) ;

template 
const double* dataPtr(const Vector<double>& vec) ;

template class LCN<double, 1>;
template class LCN<double, 2>;
template class LCN<double, 3>;
template class LCN<double, 4>;

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const Vector<double>& x);

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const LCN<double, 1>& x);

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const LCN<double, 2>& x);

template Vector<double> operator*(
  const LinearOperator<double>& A,
  const LCN<double, 3>& x);

template Vector<double>& Vector<double>::operator+=(const LCN<double, 3>& x);
template Vector<double>& Vector<double>::operator-=(const LCN<double, 3>& x);

template double norm1(const LCN<double, 1>& x);
template double norm1(const LCN<double, 2>& x);
template double norm1(const LCN<double, 3>& x);

template double norm2(const LCN<double, 1>& x);
template double norm2(const LCN<double, 2>& x);
template double norm2(const LCN<double, 3>& x);

template double normInf(const LCN<double, 1>& x);
template double normInf(const LCN<double, 2>& x);
template double normInf(const LCN<double, 3>& x);

template double min(const LCN<double, 1>& x);
template double min(const LCN<double, 2>& x);
template double min(const LCN<double, 3>& x);

template double max(const LCN<double, 1>& x);
template double max(const LCN<double, 2>& x);
template double max(const LCN<double, 3>& x);

template Vector<double> abs(const LCN<double, 1>& x);
template Vector<double> abs(const LCN<double, 2>& x);
template Vector<double> abs(const LCN<double, 3>& x);

template Vector<double> reciprocal(const LCN<double, 1>& x);
template Vector<double> reciprocal(const LCN<double, 2>& x);
template Vector<double> reciprocal(const LCN<double, 3>& x);

template LCN<double, 1> operator*(const double& a, const Vector<double>& x);
template LCN<double, 1> operator*(const Vector<double>& x, const double& a);
template LCN<double, 1> operator/(const Vector<double>& x, const double& a);

template LCN<double, 1> operator*(const double& a, const LCN<double, 1>& x);
template LCN<double, 1> operator*(const LCN<double, 1>& x, const double& a);
template LCN<double, 1> operator/(const LCN<double, 1>& x, const double& a);

template LCN<double, 2> operator*(const double& a, const LCN<double, 2>& x);
template LCN<double, 2> operator*(const LCN<double, 2>& x, const double& a);
template LCN<double, 2> operator/(const LCN<double, 2>& x, const double& a);

template LCN<double, 3> operator*(const double& a, const LCN<double, 3>& x);
template LCN<double, 3> operator*(const LCN<double, 3>& x, const double& a);
template LCN<double, 3> operator/(const LCN<double, 3>& x, const double& a);

template LCN<double, 2> 
operator+(const Vector<double>& y, const LCN<double, 1>& x);
template LCN<double, 2> 
operator+(const LCN<double, 1>& x, const Vector<double>& y);
template LCN<double, 2> 
operator+(const LCN<double, 1>& x, const LCN<double, 1>& y);

template LCN<double, 2> 
operator-(const Vector<double>& y, const LCN<double, 1>& x);
template LCN<double, 2> 
operator-(const LCN<double, 1>& x, const Vector<double>& y);
template LCN<double, 2> 
operator-(const LCN<double, 1>& x, const LCN<double, 1>& y);


template LCN<double, 3> 
operator+(const Vector<double>& y, const LCN<double, 2>& x);
template LCN<double, 3> 
operator+(const LCN<double, 2>& x, const Vector<double>& y);
template LCN<double, 3> 
operator+(const LCN<double, 2>& x, const LCN<double, 1>& y);
template LCN<double, 3> 
operator+(const LCN<double, 1>& x, const LCN<double, 2>& y);

template LCN<double, 3> 
operator-(const Vector<double>& y, const LCN<double, 2>& x);
template LCN<double, 3> 
operator-(const LCN<double, 2>& x, const Vector<double>& y);
template LCN<double, 3> 
operator-(const LCN<double, 2>& x, const LCN<double, 1>& y);
template LCN<double, 3> 
operator-(const LCN<double, 1>& x, const LCN<double, 2>& y);



}

#endif
