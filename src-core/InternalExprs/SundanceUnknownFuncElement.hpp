/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNKNOWNFUNCELEMENT_H
#define SUNDANCE_UNKNOWNFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFuncElement.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    class UnknownFunctionStub;
  }

  namespace Internal
  {
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * UnknownFuncElement represents a scalar-valued element of a (possibly)
     * list-valued UnknownFunction
     */
    class UnknownFuncElement : public SymbolicFuncElement
    {
    public:
      /** */
      UnknownFuncElement(const UnknownFunctionStub* master,
                         const string& name,
                         const string& suffix,
                         int myIndex);

      /** virtual destructor */
      virtual ~UnknownFuncElement() {;}

      /** Get the master unknown function 
       * of which this object is an element */
      const UnknownFunctionStub* master() const {return master_;}


      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      const UnknownFunctionStub* master_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
