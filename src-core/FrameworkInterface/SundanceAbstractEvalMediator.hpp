/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALMEDIATOR_H
#define SUNDANCE_EVALMEDIATOR_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceEvalVector.hpp"
#include "TSFObjectWithVerbosity.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  class CoordExpr;
  class CellDiameterExpr;

  
  namespace Internal 
  {
    class MultiIndex; 
    class DiscreteFuncElement;
  }

  using namespace Internal;

  namespace Internal
    {
      /**
       * Base class for evaluation mediator objects. 
       * Evaluation mediators are responsible
       * for evaluating those expressions whose
       * calculation must be delegated to the framework.
       */
      class AbstractEvalMediator 
        : public TSFExtended::ObjectWithVerbosity<AbstractEvalMediator>
        {
        public:
          /** */
          AbstractEvalMediator();

          /** */
          virtual ~AbstractEvalMediator(){;}

          

          /** Evaluate the given coordinate expression, putting
           * its numerical values in the given EvalVector. */
          virtual void evalCoordExpr(const CoordExpr* expr,
                                     RefCountPtr<EvalVector>& vec) const = 0 ;

          /** Evaluate the given discrete function, putting
           * its numerical values in the given EvalVector. */
          virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                               const Array<MultiIndex>& mi,
                                               Array<RefCountPtr<EvalVector> >& vec) const = 0 ;

          /** Evaluate the given cell diameter expression, putting
           * its numerical values in the given EvalVector. */
          virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                            RefCountPtr<EvalVector>& vec) const = 0 ;
        private:
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
