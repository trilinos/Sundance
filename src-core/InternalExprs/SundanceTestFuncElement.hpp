/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TESTFUNCELEMENT_H
#define SUNDANCE_TESTFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFuncElement.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace FrameworkInterface
  {
    class TestFunctionBase;
  }

  namespace Internal
  {
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * TestFuncElement represents a scalar-valued element of a (possibly)
     * list-valued TestFunction
     */
    class TestFuncElement : public SymbolicFuncElement
    {
    public:
      /** */
      TestFuncElement(const TestFunctionBase* master,
                      const string& name,
                      int myIndex);

      /** virtual destructor */
      virtual ~TestFuncElement() {;}

      /** Get the master test function 
       * of which this object is an element */
      const TestFunctionBase* master() const {return master_;}

      /** Append my func ID to the set of test IDs */
      virtual void accumulateTestSet(Set<int>& testIDs) const 
      {testIDs.put(funcID());}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      const TestFunctionBase* master_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
