/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_COORDEXPREVALUATOR_H
#define SUNDANCE_COORDEXPREVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore 
{
  class CoordExpr;

  namespace Internal 
  {
    /** 
     *
     */
    class CoordExprEvaluator : public SubtypeEvaluator<CoordExpr>
    {
    public:
      /** */
      CoordExprEvaluator(const CoordExpr* expr, 
                         const EvalContext& context);

      /** */
      virtual ~CoordExprEvaluator(){;}

      /** */
      TEUCHOS_TIMER(coordEvalTimer, "coord function evaluation");

      /** */
      virtual void internalEval(const EvalManager& mgr,
                   Array<double>& constantResults,
                   Array<RefCountPtr<EvalVector> >& vectorResults) const ;
      
    private:
      
      bool doValue_;

      bool doDeriv_;

      string stringRep_;
    };
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  



#endif
