/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_GENERICLEFTPRECONDITIONER_HPP
#define PLAYA_GENERICLEFTPRECONDITIONER_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   * A one-size-fits-most left preconditioner that can be constructed by
   * accepting an operator for the left op of the preconditioner. 
   */
  template <class Scalar>
  class GenericLeftPreconditioner : public PreconditionerBase<Scalar>
  {
  public:
    /** construct with  */
    GenericLeftPreconditioner(const LinearOperator<Scalar>& left) 
    : PreconditionerBase<Scalar>(), left_(left) {;}

    /** virtual dtor */
    virtual ~GenericLeftPreconditioner(){;}

    
    /** Return the left operator */
    virtual LinearOperator<Scalar> left() const {return left_;}

    /** A call to right() results in an error for a left precond. */
    virtual LinearOperator<Scalar> right() const
    {
      TEST_FOR_EXCEPTION(true, logic_error, "right() called for a "
                         "preconditioner known to be a left precond");
      return LinearOperator<Scalar>();
    }

    /** return true because 
     * this preconditioner has a nontrivial left component. */
    virtual bool hasLeft() const {return true;}

    /** return false, because this preconditioner has
     * no nontrivial right component */
    virtual bool hasRight() const {return false;}

    /* Handleable boilerplate */
    GET_RCP(PreconditionerBase<Scalar>);

  private:
    LinearOperator<Scalar> left_;
  };

}
