#ifndef SUNDANCE_UNARYMINUSEVALUATOR_H
#define SUNDANCE_UNARYMINUSEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore 
{
  namespace Internal 
  {
    class UnaryMinus;
    
    /**
     *
     */
    class UnaryMinusEvaluator : public UnaryEvaluator<UnaryMinus>
    {
    public:
      /** */
      UnaryMinusEvaluator(const UnaryMinus* expr,
                          const EvalContext& context);

      /** */
      virtual ~UnaryMinusEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                   Array<double>& constantResults,
                   Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      TEUCHOS_TIMER(evalTimer, "unary minus evaluation");
    private:
    }; 
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
