/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTFUNCTIONBASE_H
#define SUNDANCE_TESTFUNCTIONBASE_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace FrameworkInterface
  {
    class DiscreteFunctionBase;

    using namespace Teuchos;
    using namespace Internal;

    using std::string;
    using std::ostream;

    /** 
     * TestFunctionBase is the base class for test functions. 
     * Each framework will need to implement its own subclass of
     * TestFunctionBase. 
     *
     * The interface is left very minimal so as to not place
     * any constraints on how a framework might specify the basis.
     * When a framework needs any information about the
     * test function, it will have to get it by downcasting
     * to the appropriate framework-specific subclass.
     *
     * <h4> Writing a TestFunctionBase subclass </h4>
     *
     * For purposes of interaction with the Sundance core, no 
     * additional methods are required.
     * However, most frameworks will require extensions to 
     * TestFunctionBase that can supply the framework with information
     * on the basis used by the test func. See the
     * demo and standard frameworks for information on how to do this.
     */
    class TestFunctionBase : public SymbolicFunc
    {
    public:
      /** Construct a scalar-valued test function */
      TestFunctionBase(const string& name, int nElems=1);

      /** */
      virtual ~TestFunctionBase() {;}

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
    };
  }
}


#endif
