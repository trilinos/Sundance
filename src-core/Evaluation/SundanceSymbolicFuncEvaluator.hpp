/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SYMBOLICFUNCEVALUATOR_H
#define SUNDANCE_SYMBOLICFUNCEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore 
{
  namespace Internal 
  {
    class DiscreteFuncElementEvaluator; 
    class SymbolicFuncElement;
    /** 
     *
     */
    class SymbolicFuncElementEvaluator 
      : public SubtypeEvaluator<SymbolicFuncElement>
    {
    public:
      /** */
      SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr, 
                                   const EvalContext& context);

      /** */
      virtual ~SymbolicFuncElementEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                   Array<double>& constantResults,
                   Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      TEUCHOS_TIMER(symbolicFuncEvalTimer, "symbolic function evaluation");

      /** */
      const DiscreteFuncElementEvaluator* dfEval() const {return dfEval_;}

    private:
      Array<MultiIndex> mi_;
      Array<int> spatialDerivs_;
      Array<int> ones_;
      const DiscreteFuncElement* df_;
      const DiscreteFuncElementEvaluator* dfEval_;
      Array<string> stringReps_;
    };

    

  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  
#endif
