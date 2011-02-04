/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_LINEAROPERATORBASEDECL_HPP
#define PLAYA_LINEAROPERATORBASEDECL_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaObjectWithVerbosity.hpp"

namespace Playa
{

template <class Scalar>  class VectorSpaceBase;
template <class Scalar>  class Vector;
template <class Scalar>  class VectorType;
using Teuchos::ETransp;

/** 
 * Base class for linear operators. Most operator subtypes can safely derive
 * from LinearOpWithSpaces  which provides trivial implementations of 
 * the domain() and range() methods.
 * 
 */
template <class Scalar>
class LinearOperatorBase
  : public ObjectWithVerbosity
{
public:
  /** Virtual dtor */
  ~LinearOperatorBase(){}

  /** Return the domain */
  virtual const RCP<const VectorSpaceBase<Scalar> > domain() const = 0 ;

  /** Return the range */
  virtual const RCP<const VectorSpaceBase<Scalar> > range() const = 0 ;

  /** 
   * Apply the operator. 
   * 
   * \param applyType Indicates whether to apply the operator, its transpose,
   * or its conjugate transpose. 
   * \param in The vector on which the operator is to act
   * \param out The vector into which the result of the operation 
   * is to be written. This vector should already be initialized by the
   * appropriate space.
   **/
  virtual void apply(
    Teuchos::ETransp applyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const = 0 ;
};




}


#endif
