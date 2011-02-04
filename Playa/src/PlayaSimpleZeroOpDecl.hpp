/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_SIMPLE_ZERO_OP_DECL_HPP
#define PLAYA_SIMPLE_ZERO_OP_DECL_HPP


#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{

/** */
template <class Scalar>
class SimpleZeroOp : public LinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleZeroOp(const VectorSpace<Scalar>& domain,
    const VectorSpace<Scalar>& range);

  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  /** */
  std::string description() const ;
};


/** \relates SimpleZeroOp */
template <class Scalar>
LinearOperator<Scalar> zeroOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range);


}

#endif
