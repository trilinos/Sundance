/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_USERDEFOP_H
#define SUNDANCE_USERDEFOP_H

#include "SundanceDefs.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefOpEvaluator.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceNonlinearExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"



namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

    /**
     * UserDefOp provides a hook for inserting a user-defined nonlinear
     * function into the Sundance Expr system.
     */
  class UserDefOp : public ExprWithChildren,
                    public Internal::NonlinearExpr,
                    public GenericEvaluatorFactory<UserDefOp, UserDefOpEvaluator>
    {
    public:
      /** construct with an argument and the functor defining the operation */
      UserDefOp(const Expr& arg,
                const RefCountPtr<UserDefFunctor>& op);

      /** virtual destructor */
      virtual ~UserDefOp() {;}

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                const Set<int>& allFuncIDs,
                                bool regardFuncsAsConstant) const ;


      /** Access to the operator */
      const UserDefFunctor* op() const {return op_.get();}
    private:
      RefCountPtr<UserDefFunctor> op_;

      /** */
      static Array<RefCountPtr<ScalarExpr> > getScalarArgs(const Expr& args);

    };
}

#endif
