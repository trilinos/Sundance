/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTFUNCTIONSTUB_H
#define SUNDANCE_TESTFUNCTIONSTUB_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    class DiscreteFunctionStub;

    using namespace Teuchos;
    using namespace Internal;

    using std::string;
    using std::ostream;

    /** 
     * TestFunctionStub is the base class for test functions. 
     * Each framework will need to implement its own subclass of
     * TestFunctionStub. 
     *
     * The interface is left very minimal so as to not place
     * any constraints on how a framework might specify the basis.
     * When a framework needs any information about the
     * test function, it will have to get it by downcasting
     * to the appropriate framework-specific subclass.
     *
     * <h4> Writing a TestFunctionStub subclass </h4>
     *
     * For purposes of interaction with the Sundance core, no 
     * additional methods are required.
     * However, most frameworks will require extensions to 
     * TestFunctionStub that can supply the framework with information
     * on the basis used by the test func. See the
     * demo and standard frameworks for information on how to do this.
     */
    class TestFunctionStub : public SymbolicFunc
    {
    public:
      /** Construct a scalar-valued test function */
      TestFunctionStub(const string& name, int nElems=1);

      /** */
      virtual ~TestFunctionStub() {;}

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
    };
  }
}
                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  




#endif
