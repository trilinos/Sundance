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
   *
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
      : UnknownFunctionStub(name, basis.length()), FuncWithBasis(basis)
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
