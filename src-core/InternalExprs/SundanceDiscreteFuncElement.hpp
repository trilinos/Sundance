/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETEFUNCELEMENT_H
#define SUNDANCE_DISCRETEFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceLeafExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    class DiscreteFunctionStub;
  }

  namespace Internal
  {
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * DiscreteFuncElement represents a scalar-valued element
     * of a (possibly) vector-valued discrete function. 
     *
     * DiscreteFuncElement is framework-independent. Any framework-specific
     * information should go in a subclass of DiscreteFunctionStub.
     * The DiscreteFunctionStub object can be accessed through the
     * <tt>master()</tt> method of this class.
     */
    class DiscreteFuncElement : public virtual LeafExpr,
                                public virtual FuncElementBase
    {
    public:
      /** */
      DiscreteFuncElement(DiscreteFunctionStub* master, 
                          const string& name,
                          int myIndex);

      /** virtual destructor */
      virtual ~DiscreteFuncElement() {;}

      /** Get the master discrete function 
       * of which this object is an element */
      const DiscreteFunctionStub* master() const {return master_;}

      /** Get the master discrete function 
       * of which this object is an element */
      DiscreteFunctionStub* master() {return master_;}

      /** Get my index into the master's list of elements */
      int myIndex() const {return myIndex_;}

      /**
       * Indicate whether the given functional derivative is nonzero.
       */
      virtual bool hasNonzeroDeriv(const MultipleDeriv& d) const ;


      /**
       * Find all functions and their derivatives beneath my level
       * in the tree. For discrete functions, this does nothing. 
       */
      virtual void getRoughDependencies(Set<Deriv>& /* funcs */) const {;}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:

      DiscreteFunctionStub* master_;

      int myIndex_;
      
#endif /* DOXYGEN_DEVELOPER_ONLY */
    };
  }
}

#endif
