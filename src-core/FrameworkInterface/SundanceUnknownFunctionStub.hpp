/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNKNOWNFUNCTIONSTUB_H
#define SUNDANCE_UNKNOWNFUNCTIONSTUB_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;

  namespace Internal
  {
    using std::string;
    using std::ostream;

    /** 
   * UnknownFunctionStub is the base class for unknown functions. 
   * Each framework will need to implement its own subclass of
   * UnknownFunctionStub. 
   *
   * The interface is left very minimal so as to not place
   * any constraints on how a framework might specify the basis.
   * When a framework needs any information about the
   * unknown function, it will have to get it by downcasting
   * to the appropriate framework-specific subclass.
   *
   * <h4> Writing a UnknownFunctionStub subclass </h4>
   *
   * For purposes of interaction with the Sundance core, no 
   * additional methods are required.
   * However, most frameworks will require extensions to 
   * UnknownFunctionStub that can supply the framework with information
   * on the basis used by the unknown func. See the
   * demo and standard frameworks for information on how to do this.
   */
    class UnknownFunctionStub : public SymbolicFunc
    {
    public:
      /** */
      UnknownFunctionStub(const string& name, int nElems=1);

      /** virtual destructor */
      virtual ~UnknownFunctionStub() {;}

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}


    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
