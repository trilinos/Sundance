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
        static DerivSet setupExpr(const Expr& expr, 
                                  const Expr& tests,
                                  const Expr& unks,
                                  const Expr& u0,
                                  const RegionQuadCombo& region, 
                                  const EvaluatorFactory* factory);
        /** */
        static DerivSet setupExpr(const Expr& expr, 
                                  const RegionQuadCombo& region, 
                                  const EvaluatorFactory* factory);

        /** */
        static DerivSet identifyNonzeroDerivs(const Expr& expr,
                                              const Expr& tests,
                                              const Expr& unks,
                                              const Expr& u0);

      /** */
      static Time& preprocTimer() 
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Expr preprocessing"); 
        return *rtn;
      }
      };
  }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
