/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALUATORFACTORY_H
#define SUNDANCE_EVALUATORFACTORY_H

#include "SundanceDefs.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {
      class Evaluator;
      class EvalContext;
      class EvaluatableExpr;

      using namespace Teuchos;

      /**
       *
       */
      class EvaluatorFactory 
        {
        public:
          /** */
          EvaluatorFactory(){;}

          /** */
          virtual ~EvaluatorFactory(){;}

          /** */
          virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                             const EvalContext& context) const = 0 ;

        };
      
      
      /**
       *
       */
      template <class ExprT, class EvalT>
      class GenericEvaluatorFactory : virtual public EvaluatorFactory
        {
        public:
          /** */
          GenericEvaluatorFactory(){;}

          /** */
          virtual ~GenericEvaluatorFactory(){;}

          /** */
          virtual Evaluator* 
          createEvaluator(const EvaluatableExpr* expr,
                          const EvalContext& context) const 
            {
              return new EvalT(dynamic_cast<const ExprT*>(expr), context);
            }
        };
    }
}

#endif
