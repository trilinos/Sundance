/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DISCRETEFUNCEVALUATOR_H
#define SUNDANCE_DISCRETEFUNCEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace SundanceCore 
{
  namespace Internal 
  {
    class DiscreteFuncElement;

    /** 
     *
     */
    class DiscreteFuncElementEvaluator 
      : public SubtypeEvaluator<DiscreteFuncElement>
    {
    public:
      /** */
      DiscreteFuncElementEvaluator(const DiscreteFuncElement* expr, 
                                   const EvalContext& context);

      /** */
      virtual ~DiscreteFuncElementEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                Array<double>& constantResults,
                Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      TEUCHOS_TIMER(discreteFuncEvalTimer, "discrete function evaluation");

      /** */
      int miIndex(const MultiIndex& mi) const ;

      /** */
      bool hasMultiIndex(const MultiIndex& mi) const ;

    private:
      Array<MultiIndex> mi_;

      Map<MultiIndex, int> miToIndexMap_;

      Array<string> stringReps_;
    };

  }
}



#endif
