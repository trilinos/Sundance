/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALUATORFACTORY_H
#define SUNDANCE_EVALUATORFACTORY_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      class Evaluator;
      class EvaluatableExpr;
      using namespace Teuchos;
      
      /**
       *
       */
      class EvaluatorFactory
        {
        public:
          /** */
          EvaluatorFactory();

          /** */
          virtual ~EvaluatorFactory(){;}

          /** */
          virtual Evaluator* 
          createEvaluator(const EvaluatableExpr* expr,
                          int derivSetIndex) const = 0 ;

          /** */
          static RefCountPtr<EvaluatorFactory>& defaultEvaluator();

        protected:
          /** Method to create those evaluators that can be built by
           * the base evaluator factory */
          Evaluator* commonCreate(const EvaluatableExpr* expr) const ;
        private:
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
