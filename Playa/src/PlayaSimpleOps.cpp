/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaDefs.hpp"

#ifdef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION

#include "PlayaSimpleAddedOpImpl.hpp"
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaSimpleComposedOpImpl.hpp"
#include "PlayaSimpleDiagonalOpImpl.hpp"
#include "PlayaSimpleIdentityOpImpl.hpp"
#include "PlayaSimpleScaledOpImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaSimpleZeroOpImpl.hpp"

namespace Playa
{

template class SimpleAddedOp<double>;
template class SimpleBlockOp<double>;
template class SimpleComposedOp<double>;
template class SimpleDiagonalOp<double>;
template class SimpleIdentityOp<double>;
template class SimpleScaledOp<double>;
template class SimpleTransposedOp<double>;
template class SimpleZeroOp<double>;

template 
LinearOperator<double> operator*(const LinearOperator<double>& A,
  const LinearOperator<double>& B);

template 
LinearOperator<double> operator+(const LinearOperator<double>& A,
  const LinearOperator<double>& B);

template 
LinearOperator<double> operator*(const double& s,
  const LinearOperator<double>& A);

template 
LinearOperator<double> scaledOperator(const double& s,
  const LinearOperator<double>& A);

template 
LinearOperator<double> 
composedOperator(const Array<LinearOperator<double> >& A);

template 
LinearOperator<double> transposedOperator(const LinearOperator<double>& A);

template 
LinearOperator<double> diagonalOperator(const Vector<double>& v);

template
LinearOperator<double> identityOperator(const VectorSpace<double>& space);

template
LinearOperator<double> zeroOperator(const VectorSpace<double>& domain,
  const VectorSpace<double>& range);

template 
LinearOperator<double> makeBlockOperator(
  const VectorSpace<double>& domain,
  const VectorSpace<double>& range
  );

}


#endif
