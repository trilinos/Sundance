/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLDIAMETEREXPREVALUATOR_H
#define SUNDANCE_CELLDIAMETEREXPREVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore 
{
  class CellDiameterExpr;

  namespace Internal 
  {
    /** 
     *
     */
    class CellDiameterExprEvaluator : public SubtypeEvaluator<CellDiameterExpr>
    {
    public:
      /** */
      CellDiameterExprEvaluator(const CellDiameterExpr* expr, 
                                const EvalContext& context);

      /** */
      virtual ~CellDiameterExprEvaluator(){;}

      /** */
      TEUCHOS_TIMER(cellDiameterEvalTimer, "cell diameter evaluation");

      /** */
      virtual void internalEval(const EvalManager& mgr,
                   Array<double>& constantResults,
                   Array<RefCountPtr<EvalVector> >& vectorResults) const ;
      
    private:
      
      string stringRep_;
    };
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  




#endif
