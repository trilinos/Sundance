/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALMANAGER_H
#define SUNDANCE_EVALMANAGER_H

#include "SundanceDefs.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceAbstractEvalMediator.hpp"
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

      using namespace Internal;

      /**
       * EvalManager provides methods for interfacing to the framework
       * through an AbstractEvalMediator and managing temporary variables
       * through a TempStack.
       *
       * If no mediator is set, string evaluations will be done 
       */
      class EvalManager : public Noncopyable
        {
        public:
          /** Empty ctor */
          EvalManager();

          /** */
          void evalCoordExpr(const CoordExpr* expr,
                             RefCountPtr<EvalVector> const & result) const ;

          /** */
          void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                       const MultiIndex& mi,
                                       RefCountPtr<EvalVector> const & result) const ;

          void setMediator(const RefCountPtr<AbstractEvalMediator>& med) 
          {mediator_ = med;}

          /** Return a pointer to the mediator. We'll need the
           * mediator for computing framework-specific functions.
           */
          const AbstractEvalMediator* mediator() const {return mediator_.get();}

          /** */
          void setRegion(const EvalContext& region)
            {region_ = region;}

          /** */
          const EvalContext& getRegion() const {return region_;}

          /** */
          TempStack& stack() const {return stack_;}

        private:

          bool numericalEval() const {mediator_.get() != 0;}

          EvalContext region_;

          RefCountPtr<AbstractEvalMediator> mediator_;

          mutable TempStack stack_;
      };

    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
