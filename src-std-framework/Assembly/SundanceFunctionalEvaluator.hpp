/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCTIONALEVALUATOR_H
#define SUNDANCE_FUNCTIONALEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvaluatableExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * 
     */
    class FunctionalEvaluator 
      : public TSFExtended::ObjectWithVerbosity<FunctionalEvaluator>
    {
    public:
      /** */
      FunctionalEvaluator(const Mesh& mesh, 
                          const CellFilter& filter,
                          const Expr& expr,
                          const QuadratureFamily& quad,
                          const VerbositySetting& verb = classVerbosity());

      /** */
      double evaluate() const ;
      
      /** */
      static int& workSetSize() 
      {static int rtn = defaultWorkSetSize(); return rtn;}
      
      
    private:

      /** */
      static int defaultWorkSetSize() {return 100;}


      /** */
      Array<RefCountPtr<QuadratureEvalMediator> > mediator_;

      /** */
      const EvaluatableExpr* evalExpr_;

      /** */
      RefCountPtr<EvalManager> evalMgr_;
      
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
