/*@HEADER
//   
//@HEADER
*/


#ifndef PLAYA_PARAMETERLIST_PRECONDITIONERFACTORY_HPP
#define PLAYA_PARAMETERLIST_PRECONDITIONERFACTORY_HPP

#include "PlayaDefs.hpp"
#include "PlayaPreconditionerFactoryBase.hpp"
#include "PlayaILUKPreconditionerFactory.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaILUFactorizableOp.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaMLOperator.hpp"

namespace Playa
{
using namespace Teuchos;

/**
 * 
 */
class ParameterListPreconditionerFactory
  : public PreconditionerFactoryBase<double>
{
public:

  /** Construct with a parameter list */
  ParameterListPreconditionerFactory(const ParameterList& params)
    : PreconditionerFactoryBase<double>(), params_(params)
    {
      const std::string& pName = params_.name();
      TEST_FOR_EXCEPTION(pName != "Preconditioner", std::runtime_error,
        "expected tag=Preconditioner in parameter list " << std::endl 
        << params_);
    }

  /** virtual dtor */
  virtual ~ParameterListPreconditionerFactory(){;}
    
  /** */
  Preconditioner<double>
  createPreconditioner(const LinearOperator<double>& A) const ;

  /* Handleable boilerplate */
  GET_RCP(PreconditionerFactoryBase<double>);
private:

  ParameterList params_;
};


}

#endif




