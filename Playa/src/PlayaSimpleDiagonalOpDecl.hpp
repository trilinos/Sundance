/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_SIMPLE_DIAGONAL_OP_DECL_HPP
#define PLAYA_SIMPLE_DIAGONAL_OP_DECL_HPP



#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;





/** */
template <class Scalar>
class SimpleDiagonalOp : public LinearOpWithSpaces<Scalar>,
  public Printable
{
public:
  /** */
  SimpleDiagonalOp(const Vector<Scalar>& diag);


  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const ;

  /** */
  std::string description() const ;

  /** */
  void print(std::ostream& os) const ;

private:
  Vector<Scalar> diag_;
};

/** \relates SimpleDiagonalOp */
template <class Scalar>
LinearOperator<Scalar> diagonalOperator(
  const Vector<Scalar>& diag);

}

#endif
