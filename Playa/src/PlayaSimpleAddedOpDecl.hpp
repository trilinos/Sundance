/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_SIMPLE_ADDED_OP_DECL_HPP
#define PLAYA_SIMPLE_ADDED_OP_DECL_HPP



#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;



/**
 * Represent a sum of operators A_0 + A_1 + ... + A_n.
 */
template <class Scalar>
class SimpleAddedOp : public LinearOpWithSpaces<Scalar>
{
public:
  /** */
  SimpleAddedOp(const Array<LinearOperator<Scalar> >& ops);

  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;
  
  /** */
  std::string description() const ;

private:
    Array<LinearOperator<Scalar> > ops_;
};

/** \relates LinearOperator */
template <class Scalar>
LinearOperator<Scalar> addedOperator(
  const Array<LinearOperator<Scalar> >& ops);


/** \relates LinearOperator */
template <class Scalar> 
LinearOperator<Scalar> 
operator+(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B);


}

#endif
