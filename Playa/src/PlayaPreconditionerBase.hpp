/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_PRECONDITIONERBASE_HPP
#define PLAYA_PRECONDITIONERBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
  using namespace Teuchos;

  /**
   * Base class for preconditioners. A general preconditioner object
   * is split into a left preconditioner M1^-1 and a right
   * preconditioner M2^-1. To solve A x = b, we define the auxiliary
   * system M2^-1 y = x, and solve M1^-1 A M2^-1 y = M1^-1 b to obtain y.
   * Having y, we can quickly recover x by applying M2^-1 to y.
   *
   * The base class implements neither a left nor a right preconditioner.
   */
  template <class Scalar>
  class PreconditionerBase : public Playa::Handleable<PreconditionerBase<Scalar> >
  {
  public:
    /** empty ctor */
    PreconditionerBase() {;}

    /** virtual dtor */
    virtual ~PreconditionerBase(){;}

    
    /** */
    virtual LinearOperator<Scalar> left() const = 0 ;

    /** */
    virtual LinearOperator<Scalar> right() const = 0 ;

    /** return true if this preconditioner has a nontrivial left component */
    virtual bool hasLeft() const = 0 ;

    /** return true if this preconditioner has
     * a nontrivial right component */
    virtual bool hasRight() const = 0 ;

  private:
  };

}

#endif
