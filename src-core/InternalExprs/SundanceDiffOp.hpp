/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DIFFOP_H
#define SUNDANCE_DIFFOP_H

#include "SundanceDefs.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceDiffOpEvaluator.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     *
     */
    class DiffOp : public UnaryExpr,
                   public GenericEvaluatorFactory<DiffOp, DiffOpEvaluator>
    {
    public:
      /** ctor */
      DiffOp(const MultiIndex& op, const RefCountPtr<ScalarExpr>& arg);

      /** virtual destructor */
      virtual ~DiffOp() {;}

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;

      

      /** */
      const Deriv& myCoordDeriv() const {return myCoordDeriv_;}

      /** */
      const MultiIndex& mi() const {return mi_;}

      /** Get the functions that are required in the evaluation
       * of the multiple deriv d */
      const SundanceUtils::Set<Deriv>& requiredFunctions(const MultipleDeriv& d) const 
      {return requiredFunctions_[d];}

      /** */
      bool requiresFunctionsToEval(const MultipleDeriv& d) const 
      {return requiredFunctions_.containsKey(d);}

    
      
      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                const Set<int>& allFuncIDs,
                                bool regardFuncsAsConstant) const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}


      /** Determine the functional derivatives of this expression that
       * are produced through the chain rule
       * by the presence of the given derivative of the argument. 
       */
      void getResultDerivs(const MultipleDeriv& argDeriv,
                           const DerivState& sourceState,
                           Map<MultipleDeriv, DerivState>& isolatedTerms,
                           Map<MultipleDeriv, Deriv>& funcTerms) const ;


      /** */
      virtual Set<MultiIndex> argMultiIndices(const Set<MultiIndex>& myMultiindices) const ;
        
      /** */
      bool ignoreFuncTerms() const {return ignoreFuncTerms_;}



    protected:

      /** 
       * Given a set of active function combinations, get the active
       * funcs for the next higher order of differentiation.
       */
      Set<MultiSet<int> > 
      argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs,
                     const Set<int>& allFuncIDs) const ;

    private:

      
      bool canBackOutDeriv(const Deriv& d, const MultiIndex& x, 
                           Deriv& rtnDeriv) const ;

      

      MultiIndex mi_;


      Deriv myCoordDeriv_;

      mutable Map<MultipleDeriv, SundanceUtils::Set<Deriv>, 
                  increasingOrder<MultipleDeriv> > requiredFunctions_;


      mutable bool ignoreFuncTerms_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
