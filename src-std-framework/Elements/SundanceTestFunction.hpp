/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTFUNCTION_H
#define SUNDANCE_TESTFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceTestFunctionStub.hpp"
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
  class TestFunction : public TestFunctionStub,
                          public FuncWithBasis
  {
  public:
    /** */
    TestFunction(const BasisFamily& basis, const string& name="")
      : TestFunctionStub(name, basis.dim()), FuncWithBasis(basis)
    {;}

    /** */
    TestFunction(const Array<BasisFamily>& basis, const string& name="")
      : TestFunctionStub(name, BasisFamily::size(basis)), 
        FuncWithBasis(basis)
    {;}

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** virtual destructor */
    virtual ~TestFunction() {;}

    /* boilerplate */
    GET_RCP(ExprBase);
#endif /* DOXYGEN_DEVELOPER_ONLY */
  };

}



#endif
