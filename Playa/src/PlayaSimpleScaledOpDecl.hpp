/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_SIMPLE_SCALED_OP_DECL_HPP
#define PLAYA_SIMPLE_SCALED_OP_DECL_HPP



#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;




/**
 * Represent a scaled operator
 */
template <class Scalar>
class SimpleScaledOp : public LinearOpWithSpaces<Scalar>,
                       public Printable
{
public:
  /** */
  SimpleScaledOp(const Scalar& alpha, const LinearOperator<Scalar>& A);
  
  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  
  /** */
  std::string description() const ;

  /** */
  const Scalar& alpha() const {return alpha_;}

  /** */
  LinearOperator<Scalar> op() const {return A_;}

  /** */
  void print(std::ostream& os) const ;

private:
  Scalar alpha_;
  LinearOperator<Scalar> A_;
};


/** \relates SimpleScaledOp */
template <class Scalar>
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op);


template <class Scalar> 
LinearOperator<Scalar> operator*(const Scalar& a, const LinearOperator<Scalar>& A);

}

#endif
