/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_NONLINEARUNARYOPEVALUATOR_H
#define SUNDANCE_NONLINEARUNARYOPEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace SundanceCore 
{
  namespace Internal 
  {
    class NonlinearUnaryOp;
    class SymbolicFuncElementEvaluator;
    
    /**
     *
     */
    class NonlinearUnaryOpEvaluator : public UnaryEvaluator<NonlinearUnaryOp>
    {
    public:
      /** */
      NonlinearUnaryOpEvaluator(const NonlinearUnaryOp* expr,
                                const EvalContext& context);

      /** */
      virtual ~NonlinearUnaryOpEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                   Array<double>& constantResults,
                   Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      TEUCHOS_TIMER(evalTimer, "nonlinear unary op evaluation");
    private:
      int maxOrder_;
      int d0ResultIndex_;
      int d0ArgDerivIndex_;
      bool d0ArgDerivIsConstant_;
      Array<int> d1ResultIndex_;
      Array<int> d1ArgDerivIndex_;
      Array<int> d1ArgDerivIsConstant_;
      Array<int> d2ResultIndex_;
      Array<Array<int> > d2ArgDerivIndex_;
      Array<Array<int> > d2ArgDerivIsConstant_;
    }; 
  }
}


#endif
