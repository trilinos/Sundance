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
  namespace Internal
  {
    class TestFunctionStub;
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
      TestFuncElement(const TestFunctionStub* master,
                      const string& name,
                      const string& suffix,
                      int myIndex);

      /** virtual destructor */
      virtual ~TestFuncElement() {;}

      /** Get the master test function 
       * of which this object is an element */
      const TestFunctionStub* master() const {return master_;}

      /** Append my func ID to the set of test IDs */
      virtual void accumulateTestSet(Set<int>& testIDs) const 
      {testIDs.put(funcID());}



      /** Test whether all terms have test functions. 
       * I'm a test function, so return true */
      virtual bool allTermsHaveTestFunctions() const {return true;}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      const TestFunctionStub* master_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
