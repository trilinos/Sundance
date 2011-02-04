/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_PRECONDITIONERFACTORYBASE_HPP
#define PLAYA_PRECONDITIONERFACTORYBASE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaPreconditioner.hpp"
#include "Teuchos_ParameterList.hpp"


namespace Playa
{
using namespace Teuchos;

/**
 * Base class for preconditioner factories. 
 */
template <class Scalar>
class PreconditionerFactoryBase 
  : public Playa::Handleable<PreconditionerFactoryBase<Scalar> >
{
public:
  /** empty ctor */
  PreconditionerFactoryBase() {;}

  /** virtual dtor */
  virtual ~PreconditionerFactoryBase(){;}

    
  /** */
  virtual Preconditioner<Scalar> createPreconditioner(const LinearOperator<Scalar>& A) const = 0 ;

private:
};

}

#endif
