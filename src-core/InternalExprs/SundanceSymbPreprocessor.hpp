/* @HEADER@ */
/* @HEADER@ */



#ifndef SUNDANCE_SYMBPREPROCESSOR_H
#define SUNDANCE_SYMBPREPROCESSOR_H

#include "SundanceExpr.hpp"
#include "SundanceEvaluatableExpr.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;

  using std::string;
  using std::ostream;

  namespace Internal
    {
      /** */
      class SymbPreprocessor
      {
      public:

        /** */
        static DerivSet setupVariations(const Expr& expr, 
                                        const Expr& vars,
                                        const Expr& varEvalPts,
                                        const Expr& unks,
                                        const Expr& unkEvalPts,
                                        const Expr& fixedFields,
                                        const Expr& fixedFieldEvalPts,
                                        const EvalContext& region);

        /** */
        static DerivSet setupFunctional(const Expr& expr, 
                                        const Expr& fixedFields,
                                        const Expr& fixedFieldEvalPts,
                                        const EvalContext& region);

        /** */
        static DerivSet setupGradient(const Expr& expr, 
                                      const Expr& vars,
                                      const Expr& varEvalPts,
                                      const Expr& fixedFields,
                                      const Expr& fixedFieldEvalPts,
                                      const EvalContext& region);

        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const Expr& tests,
                                  const Expr& unks,
                                  const Expr& u0,
                                  const EvalContext& region);
        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const Expr& unks,
                                  const Expr& u0,
                                  const EvalContext& region);
        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const EvalContext& region);

        /** */
        TEUCHOS_TIMER(preprocTimer, "symbolic preprocessing");
      };
  }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
