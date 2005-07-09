#ifndef SUNDANCE_DIFFOPEVALUATOR_H
#define SUNDANCE_DIFFOPEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore 
{
  namespace Internal 
  {
    class DiffOp;
    class DiscreteFuncElementEvaluator;
    
    /**
     *
     */
    class DiffOpEvaluator : public UnaryEvaluator<DiffOp>
    {
    public:
      /** */
      DiffOpEvaluator(const DiffOp* expr,
                      const EvalContext& context);

      /** */
      virtual ~DiffOpEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                Array<double>& constantResults,
                Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** We need a specialized resetting method for diff op
       * evaluators that also resets the discrete func evaluators
       * used in the functional chain rule */
      virtual void resetNumCalls() const ;

      /** */
      TEUCHOS_TIMER(evalTimer, "diff op evaluation");
    private:

      Array<int> isConstant_;

      Array<int> resultIndices_;
      
      Array<Array<int> > constantMonomials_;

      Array<Array<int> > vectorMonomials_;

      Array<Array<int> > constantFuncCoeffs_;

      Array<Array<int> > vectorFuncCoeffs_;

      Array<const DiscreteFuncElementEvaluator*> funcEvaluators_;

      /** Indices into the function evaluator table for the funcs
       * appearing with constant coeffs in the chain rule */
      Array<Array<int> > constantCoeffFuncIndices_;

      /** Indices into the list of multiindices for the funcs
       * appearing with constant coeffs in the chain rule */
      Array<Array<int> > constantCoeffFuncMi_;

      /** Indices into the function evaluator table for the funcs
       * appearing with vector coeffs in the chain rule */
      Array<Array<int> > vectorCoeffFuncIndices_;

      /** Indices into the list of multiindices for the funcs
       * appearing with vector coeffs in the chain rule */
      Array<Array<int> > vectorCoeffFuncMi_;
    }; 
  }
}
                 
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
