/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SYMBOLICFUNCELEMENT_H
#define SUNDANCE_SYMBOLICFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    class DiscreteFuncElement;
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * SymbolicFuncElement represents a scalar-valued element of a (possibly)
     * list-valued SymbolicFunction. 
     */
    class SymbolicFuncElement : public FuncElementBase,
                                virtual public EvaluatableExpr,
                                public GenericEvaluatorFactory<SymbolicFuncElement, SymbolicFuncElementEvaluator>
    {
    public:
      /** */
      SymbolicFuncElement(const string& name, const string& suffix,
                          int myIndex);
      
      /** virtual destructor */
      virtual ~SymbolicFuncElement() {;}

      /** Get my index into the master's list of elements */
      int myIndex() const {return myIndex_;}


      /** Specify that expressions involving this function are to be evaluated
       * with this function set to zero. Test functions should always be
       * evaluated at zero. For unknown functions, 
       * substituting zero is appropriate for computing
       * the functional derivatives that arise in a linear problem.
       * */
      void substituteZero() const ;

      /** Specify that expressions involving this function are to be evaluated
       * with this function set to the discrete function \f$u_0\f$. 
       * This is appropriate for computing
       * the functional derivatives that arise in a nonlinear expression
       * being linearized about \f$u_0\f$. 
       */
      void substituteFunction(const RefCountPtr<DiscreteFuncElement>& u0) const ;

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      const EvaluatableExpr* evalPt() const {return evalPt_.get();}

      

      /** \name Preprocessing */
      //@{
      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which functional derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const ;

      

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      mutable RefCountPtr<EvaluatableExpr> evalPt_;

      mutable Array<int> evalPtDerivSetIndices_;

      int myIndex_;
      

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
