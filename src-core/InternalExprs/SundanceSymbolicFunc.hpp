/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SYMBOLICFUNC_H
#define SUNDANCE_SYMBOLICFUNC_H


#include "SundanceDefs.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {

    using namespace Internal;
    using namespace Teuchos;
    
    using std::string;
    using std::ostream;

    /** 
     * SymbolicFunc is a base class for functions such as test and unknown
     * functions that are "variables" in a weak form. Symbolic functions
     * cannot be evaluated directly; before evaluating a weak form,
     a value must be substituted for 
     * each symbolic func using either the substituteZero() or
     * substituteFunction() method. 
     *
     * A symbolic function will construct itself as a list of
     * SymbolicFuncElement objects that point back to the SymbolicFunction.
     */
    class SymbolicFunc : public ListExpr
    {
    public:
      /** Empty ctor, initializes list to empty */
      SymbolicFunc();

      /** virtual destructor */
      virtual ~SymbolicFunc() {;}

      /** Specify that expressions involving this function are to be evaluated
       * with this function set to zero. This is appropriate for computing
       * the functional derivatives that arise in a linear problem */
      void substituteZero() const ;

      /** Specify that expressions involving this function are to be evaluated
       * with this function set to the discrete function \f$u_0\f$. 
       * This is appropriate for computing
       * the functional derivatives that arise in a nonlinear expression
       * being linearized about \f$u_0\f$. 
       */
      void substituteFunction(const RefCountPtr<DiscreteFunctionStub>& u0) const ;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
