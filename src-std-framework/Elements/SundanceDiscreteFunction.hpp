/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETEFUNCTION_H
#define SUNDANCE_DISCRETEFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceFuncWithBasis.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Internal;

  /** 
   *
   */
  class DiscreteFunction : public DiscreteFunctionStub,
                          public FuncWithBasis
  {
  public:
    /** */
    DiscreteFunction(const BasisFamily& basis, const string& name="")
      : DiscreteFunctionStub(name, basis.dim()), FuncWithBasis(basis)
    {;}

    /** */
    DiscreteFunction(const Array<BasisFamily>& basis, const string& name="")
      : DiscreteFunctionStub(name, BasisFamily::size(basis)), 
        FuncWithBasis(basis)
    {;}

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** virtual destructor */
    virtual ~DiscreteFunction() {;}

    /* boilerplate */
    GET_RCP(ExprBase);
#endif /* DOXYGEN_DEVELOPER_ONLY */
  };

}



#endif
