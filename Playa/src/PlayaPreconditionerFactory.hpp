/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_PRECONDITIONERFACTORY_HPP
#define PLAYA_PRECONDITIONERFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaPreconditionerFactoryBase.hpp"

namespace Playa
{
  /**
   * PreconditionerFactory builds an implementation-specific preconditioner
   * from an abstract specification. 
   * 
   * Preconditioners are constructed indirectly through factories
   * rather then directly by preconditioner ctor calls. 
   * The reason for this is that when we create a solver and want to 
   * specify the preconditioner, we don't yet know the matrix (or even
   * the type or size of matrix) on which the solver is going to operate.
   * Thus we have to defer construction of the preconditioner until 
   * the solve() call when the matrix is available. The factory gives
   * us a means by which we can build a preconditioner at that point.
   */
  template <class Scalar> 
  class PreconditionerFactory 
    : public Playa::Handle<PreconditionerFactoryBase<Scalar> >
  {
  public:
    /* Boilerplate ctors */
    HANDLE_CTORS(PreconditionerFactory<Scalar>, PreconditionerFactoryBase<Scalar>);

    /** create a concrete preconditioner */
    Preconditioner<Scalar> createPreconditioner(const LinearOperator<Scalar>& A) const
    {return this->ptr()->createPreconditioner(A);}
    
  };
}

#endif
