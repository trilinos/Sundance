/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_SIMPLE_IDENTITY_OP_DECL_HPP
#define PLAYA_SIMPLE_IDENTITY_OP_DECL_HPP



#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;





/** */
template <class Scalar>
class SimpleIdentityOp : public LinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleIdentityOp(const VectorSpace<Scalar>& space);


  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const ;

  /** */
  std::string description() const ;
};

/** \relates SimpleIdentityOp */
template <class Scalar>
LinearOperator<Scalar> identityOperator(
  const VectorSpace<Scalar>& space);

}

#endif
