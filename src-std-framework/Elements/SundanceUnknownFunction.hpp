/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNKNOWNFUNCTION_H
#define SUNDANCE_UNKNOWNFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceFuncWithBasis.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Internal;

  /** 
   * UnknownFunction represents an unknown function in a finite
   * element problem. Unknown functions can be scalar or vector valued, as
   * determined at runtime through the type of basis with which
   * they are constructed.
   */
  class UnknownFunction : public UnknownFunctionStub,
                          public FuncWithBasis
  {
  public:
    /** */
    UnknownFunction(const BasisFamily& basis, const string& name="")
      : UnknownFunctionStub(name, basis.dim()), FuncWithBasis(basis)
    {;}

    /** */
    UnknownFunction(const Array<BasisFamily>& basis, const string& name="")
      : UnknownFunctionStub(name, BasisFamily::size(basis)), 
        FuncWithBasis(basis)
    {;}

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** virtual destructor */
    virtual ~UnknownFunction() {;}

    /* boilerplate */
    GET_RCP(ExprBase);
#endif /* DOXYGEN_DEVELOPER_ONLY */
  };

}



#endif
