/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALMANAGER_H
#define SUNDANCE_EVALMANAGER_H

#include "SundanceDefs.hpp"
#include "SundanceEvalRegion.hpp"
#include "SundanceEvalMediator.hpp"
#include "SundanceTempStack.hpp"
#include "SundanceNoncopyable.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  class CoordExpr;


  namespace Internal
    {
      class MultiIndex;
      class DiscreteFuncElement;

      using namespace FrameworkInterface;

      /**
       * EvalManager provides methods for interfacing to the framework
       * through an EvalMediator and managing temporary variables
       * through a TempStack.
       */
      class EvalManager : public Noncopyable
        {
        public:
          /** Construct with an EvalMediator that knows how to create
           * vectors and evaluate framework-specific functions. */
          EvalManager(const RefCountPtr<EvalMediator>& mediator);

          /** Empty ctor, creates a null mediator and does string
           * evaluations */
          EvalManager();

          /** */
          void evalCoordExpr(const CoordExpr* expr,
                             RefCountPtr<EvalVector> const & result) const ;

          /** */
          void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                       const MultiIndex& mi,
                                       RefCountPtr<EvalVector> const & result) const ;

          /** Return a pointer to the mediator. We'll need the
           * mediator for computing framework-specific functions.
           */
          const EvalMediator* mediator() const {return mediator_.get();}

          /** */
          void setRegion(const EvalRegion& region)
            {region_ = region;}

          /** */
          const EvalRegion& getRegion() const {return region_;}

          /** */
          TempStack& stack() const {return stack_;}

        private:

          bool numericalEval_;

          EvalRegion region_;

          RefCountPtr<EvalMediator> mediator_;

          mutable TempStack stack_;
      };

    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
