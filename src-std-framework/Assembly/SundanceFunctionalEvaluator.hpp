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

  /** */
  double evaluateIntegral(const Mesh& mesh, const Expr& expr);

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
                          const Expr& integral);
      /** */
      FunctionalEvaluator(const Mesh& mesh, 
                          const Expr& integral,
                          const Expr& bcs,
                          const Expr& fields,
                          const Expr& fieldValues);
      /** */
      FunctionalEvaluator(const Mesh& mesh, 
                          const Expr& integral,
                          const Expr& bcs,
                          const Expr& vars,
                          const Expr& varEvalPts,
                          const Expr& fields,
                          const Expr& fieldValues);


      /** */
      double evaluate() const ;
      
    private:
      
    /** */
    RefCountPtr<Assembler> assembler_;
      
    };
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
