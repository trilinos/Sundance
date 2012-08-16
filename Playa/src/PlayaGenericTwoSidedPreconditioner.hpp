/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_GENERICTWOSIDEDPRECONDITIONER_HPP
#define PLAYA_GENERICTWOSIDEDPRECONDITIONER_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaPreconditionerBase.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   * A one-size-fits-most left preconditioner that can be constructed by
   * accepting an operator for the left op of the preconditioner. 
   */
  template <class Scalar>
  class GenericTwoSidedPreconditioner : public PreconditionerBase<Scalar>
  {
  public:
    /** construct with an operator for the two-sided preconditioner */
    GenericTwoSidedPreconditioner(const LinearOperator<Scalar>& left, const LinearOperator<Scalar>& right) 
    : PreconditionerBase<Scalar>(), left_(left), right_(right) {;}

    /** virtual dtor */
    virtual ~GenericTwoSidedPreconditioner(){;}

    
    /** Return the left operator */
    virtual LinearOperator<Scalar> left() const {return left_;}

    /** Return the right operator */
    virtual LinearOperator<Scalar> right() const {return right_;}

    /** return true because 
     * this preconditioner has a nontrivial left component. */
    virtual bool hasLeft() const {return true;}

    /** return true because 
     * this preconditioner has a nontrivial right component. */
    virtual bool hasRight() const {return true;}

    /* Handleable boilerplate */
    GET_RCP(PreconditionerBase<Scalar>);

  private:
    LinearOperator<Scalar> left_;
    LinearOperator<Scalar> right_;
  };

}

#endif
