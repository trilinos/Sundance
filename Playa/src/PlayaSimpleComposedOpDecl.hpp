/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_SIMPLE_COMPOSED_OP_DECL_HPP
#define PLAYA_SIMPLE_COMPOSED_OP_DECL_HPP


#include "PlayaDefs.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"


namespace Playa
{
using namespace Teuchos;




/**
 * Represent a composed operator A_0 * A_1 * ... * A_n.
 */
template <class Scalar>
class SimpleComposedOp : public LinearOpWithSpaces<Scalar>,
                         public Printable
{
public:
  /** */
  SimpleComposedOp(const Array<LinearOperator<Scalar> >& ops);
  
  /** */
  void apply(Teuchos::ETransp transApplyType,
    const Vector<Scalar>& in,
    Vector<Scalar> out) const;

  
  /** */
  std::string description() const ;

  /** */
  void print(std::ostream& os) const ;

private:
  Array<LinearOperator<Scalar> > ops_;
};


template <class Scalar>
LinearOperator<Scalar> composedOperator(
  const Array<LinearOperator<Scalar> >& ops);


template <class Scalar> 
LinearOperator<Scalar> operator*(const LinearOperator<Scalar>& A, 
  const LinearOperator<Scalar>& B);

}

#endif
