/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_NONLINEARUNARYOP_H
#define SUNDANCE_NONLINEARUNARYOP_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceNonlinearUnaryOpEvaluator.hpp"
#include "SundanceNonlinearExpr.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     *
     */
    class NonlinearUnaryOp : public UnaryExpr,
                             public NonlinearExpr,
                             public GenericEvaluatorFactory<NonlinearUnaryOp, NonlinearUnaryOpEvaluator>
    {
    public:
      /** construct with an argument and the functor defining the operation */
      NonlinearUnaryOp(const RefCountPtr<ScalarExpr>& arg, 
                       const RefCountPtr<UnaryFunctor>& op);

      /** virtual destructor */
      virtual ~NonlinearUnaryOp() {;}

      /** Preprocessing step to determine which functional 
       * derivatives are nonzero */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const ;

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}


      /** Access to the operator */
      const UnaryFunctor* op() const {return op_.get();}

      /** 
       * Given a set of active function combinations, get the active
       * function combinations requested of the argument
       */
      virtual Set<MultiSet<int> > 
      argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs) const ;
    private:

      RefCountPtr<UnaryFunctor> op_;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
