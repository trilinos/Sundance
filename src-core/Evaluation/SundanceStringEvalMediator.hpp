/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_STRINGEVALMEDIATOR_H
#define SUNDANCE_STRINGEVALMEDIATOR_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceAbstractEvalMediator.hpp"


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
       *
       */
      class StringEvalMediator : public AbstractEvalMediator
        {
        public:
          /** */
          StringEvalMediator();

          /** */
          virtual ~StringEvalMediator(){;}

          /** Evaluate the given coordinate expression, putting
           * its numerical values in the given LoadableVector. */
          virtual void evalCoordExpr(const CoordExpr* expr,
                                     RefCountPtr<EvalVector>& vec) const ;

          /** Evaluate the given discrete function, putting
           * its numerical values in the given LoadableVector. */
          virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                               const Array<MultiIndex>& mi,
                                               Array<RefCountPtr<EvalVector> >& vec) const ;

          /** Evaluate the given cell diameter expression, putting
           * its numerical values in the given EvalVector. */
          virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
                                            RefCountPtr<EvalVector>& vec) const ;
            

        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
