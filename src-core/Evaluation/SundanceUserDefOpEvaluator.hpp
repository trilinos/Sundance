/* @HEADER@ */
/* @HEADER@ */


#ifndef SUNDANCE_USERDEFOPEVALUATOR_H
#define SUNDANCE_USERDEFOPEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace SundanceCore 
{
  class UserDefOp;

  namespace Internal 
  {

    class SymbolicFuncElementEvaluator;
    
    /**
     *
     */
    class UserDefOpEvaluator : public SubtypeEvaluator<UserDefOp>
    {
    public:
      /** */
      UserDefOpEvaluator(const UserDefOp* expr,
                         const EvalContext& context);

      /** */
      virtual ~UserDefOpEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                                Array<double>& constantResults,
                                Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      virtual void resetNumCalls() const ;

      /** */
      TEUCHOS_TIMER(evalTimer, "user defined nonlinear op evaluation");

    protected:
      const RefCountPtr<SparsitySuperset>& childSparsity(int i) const
      {return childSparsity_[i];}

      const EvaluatableExpr* childExpr(int i) const
      {return childExpr_[i];}

      const RefCountPtr<Evaluator>& childEval(int i) const
      {return childEval_[i];}

      void evalChildren(const EvalManager& mgr,
                        Array<Array<double> >& constResults,
                        Array<Array<RefCountPtr<EvalVector> > >& vecResults) const ;

      void evalOperator(int numPoints,
                        const Array<double>& constantArg,
                        const Array<RefCountPtr<EvalVector> >& vectorArg,
                        const Array<int>& constantArgPtr,
                        const Array<int>& vectorArgPtr,
                        RefCountPtr<EvalVector>& opResults) const ;

      
    private:
      Array<const EvaluatableExpr*> childExpr_;

      Array<RefCountPtr<SparsitySuperset> > childSparsity_;

      Array<RefCountPtr<Evaluator> > childEval_;

      int maxOrder_;
      int d0ResultIndex_;
      Array<int> d0ArgDerivIndex_;
      Array<int> d0ArgDerivIsConstant_;
      Array<int> constantArgPtr_;
      Array<int> vectorArgPtr_;
    }; 
  }
}


#endif
