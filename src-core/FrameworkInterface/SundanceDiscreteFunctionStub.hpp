/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETEFUNCTIONSTUB_H
#define SUNDANCE_DISCRETEFUNCTIONSTUB_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceListExpr.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    using namespace Internal;
    using std::string;
    using std::ostream;

    /** 
     * DiscreteFunctionStub is the base class for discrete functions. 
     * Each framework will need to implement its own subclass of
     * DiscreteFunctionStub. 
     *
     * The interface is left very minimal so as to not place
     * any constraints on how a framework might specify vectors
     * and bases. When a framework needs any information about the
     * discrete function, it will have to get it by downcasting
     * to the appropriate framework-specific subclass.
     *
     * <h4> Writing a DiscreteFunctionStub subclass </h4>
     *
     * For purposes of interaction with the Sundance core, no 
     * additional methods are required.
     * However, most frameworks will require extensions to 
     * DiscreteFunctionStub that can supply the framework with information
     * on the basis and vector used by the discrete func. See the
     * demo and standard frameworks for information on how to do this.
     */
    class DiscreteFunctionStub : public ListExpr
    {
    public:
      /** */
      DiscreteFunctionStub(const string& name, int nElems=1);

      /** virtual destructor */
      virtual ~DiscreteFunctionStub() {;}

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}

    protected:
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
