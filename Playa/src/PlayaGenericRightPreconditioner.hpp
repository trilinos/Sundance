/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_GENERICRIGHTPRECONDITIONER_HPP
#define PLAYA_GENERICRIGHTPRECONDITIONER_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaPreconditionerBase.hpp"


namespace Playa
{
using namespace Teuchos;

/**
 * A one-size-fits-most right preconditioner that can be constructed by
 * accepting an operator for the right op of the preconditioner. 
 */
template <class Scalar>
class GenericRightPreconditioner : public PreconditionerBase<Scalar>
{
public:
  /** construct with an operator for the right preconditioner */
  GenericRightPreconditioner(const LinearOperator<Scalar>& right) 
    : PreconditionerBase<Scalar>(), right_(right) {;}

  /** virtual dtor */
  virtual ~GenericRightPreconditioner(){;}

    
  /** Return the right operator */
  virtual LinearOperator<Scalar> right() const {return right_;}

  /** A call to left() results in an error for a right precond. */
  virtual LinearOperator<Scalar> left() const
    {
      TEST_FOR_EXCEPTION(true, std::logic_error, "left() called for a "
        "preconditioner known to be a right precond");
      return LinearOperator<Scalar>();
    }

  /** return true because 
   * this preconditioner has a nontrivial right component. */
  virtual bool hasRight() const {return true;}

  /** return false, because this preconditioner has
   * no nontrivial left component */
  virtual bool hasLeft() const {return false;}

  /* Handleable boilerplate */
  GET_RCP(PreconditionerBase<Scalar>);

private:
  LinearOperator<Scalar> right_;
};



}

#endif
