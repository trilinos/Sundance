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
  namespace FrameworkInterface
  {
    class UnknownFunctionBase;
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
      UnknownFuncElement(const UnknownFunctionBase* master,
                         const string& name,
                         int myIndex);

      /** virtual destructor */
      virtual ~UnknownFuncElement() {;}

      /** Get the master unknown function 
       * of which this object is an element */
      const UnknownFunctionBase* master() const {return master_;}

      /** Append my func ID to the set of unk IDs */
      virtual void accumulateUnkSet(Set<int>& unkIDs) const 
      {unkIDs.put(funcID());}


      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      const UnknownFunctionBase* master_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
