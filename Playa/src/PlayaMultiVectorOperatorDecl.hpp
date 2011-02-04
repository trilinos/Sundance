/* @HEADER@ */

 /* @HEADER@ */

#ifndef Playa_MULTI_VECTOR_OPERATOR_DECL_HPP
#define Playa_MULTI_VECTOR_OPERATOR_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaRowAccessibleOp.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "PlayaVectorDecl.hpp"

namespace Playa
{
/** 
 * A MultiVectorOperator is a linear operator whose
 * rows or columns are represented as a multivector 
 */
template <class Scalar> 
class MultiVectorOperator 
  : public LinearOpWithSpaces<Scalar>,
    public RowAccessibleOp<Scalar>
{
public:

  /**
   * Construct from an array of vectors and a specifier for the 
   * domain space. 
   */
  MultiVectorOperator(const Teuchos::Array<Vector<Scalar> >& cols,
    const VectorSpace<Scalar>& domain);

  /** Virtual dtor */
  virtual ~MultiVectorOperator(){;}


  /** 
   * Apply does an element-by-element multiply between the input 
   * vector, x, and the diagonal values.
   */
  virtual void apply(
    Teuchos::ETransp transType,
    const Vector<Scalar>& in, 
    Vector<Scalar> out) const ;


  /** Return the kth row  */
  void getRow(const int& k, 
    Teuchos::Array<int>& indices, 
    Teuchos::Array<Scalar>& values) const ;

private:

  Teuchos::Array<Vector<Scalar> > cols_;
};

/** \relates MultiVectorOperator */
template <class Scalar> 
LinearOperator<Scalar> multiVectorOperator(
  const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain);


}

#endif
