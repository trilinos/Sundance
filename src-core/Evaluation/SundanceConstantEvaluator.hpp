/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CONSTANTEVALUATOR_H
#define SUNDANCE_CONSTANTEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore 
{
  namespace Internal 
  {
    class SpatiallyConstantExpr;

    /** 
     *
     */
    class ConstantEvaluator : public SubtypeEvaluator<SpatiallyConstantExpr>
    {
    public:
      /** */
      ConstantEvaluator(const SpatiallyConstantExpr* expr, 
                        const EvalContext& context);

      /** */
      virtual ~ConstantEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                                Array<double>& constantResults,
                                Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      TEUCHOS_TIMER(constantEvalTimer, "constant expr evaluation");
      
    }; 
  }
}
               
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
