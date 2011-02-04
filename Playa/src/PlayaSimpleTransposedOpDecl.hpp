/* @HEADER@ */
//
 /* @HEADER@ */

#ifndef PLAYA_SIMPLE_TRANSPOSED_OP_DECL_HPP
#define PLAYA_SIMPLE_TRANSPOSED_OP_DECL_HPP



#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;




/**
 * Represent the transpose of an operator
 */
template <class Scalar>
class SimpleTransposedOp : public LinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleTransposedOp(const LinearOperator<Scalar>& A);
  
  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  
  /** */
  std::string description() const ;

  /** */
  LinearOperator<Scalar> op() const {return A_;}

private:
  LinearOperator<Scalar> A_;
};


/** \relates SimpleTransposedOp */
template <class Scalar>
LinearOperator<Scalar> transposedOperator(
  const LinearOperator<Scalar>& op);

}

#endif
