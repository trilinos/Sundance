/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNARYEXPR_H
#define SUNDANCE_UNARYEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExprWithChildren.hpp"
#include "SundanceExpr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using std::string;

  namespace Internal
    {
      /** 
       * UnaryExpr is a base class for operators that map a single
       * scalar expr to another scalar expr. 
       */
      class UnaryExpr : public ExprWithChildren
        {
        public:
          /** construct with the argument */
          UnaryExpr(const RefCountPtr<ScalarExpr>& arg);

          /** virtual dtor */
          virtual ~UnaryExpr() {;}

          /** Return the operand */
          Expr arg() const {return child(0);}

          /** Downcast the argument to an evaluatable expr */
          const EvaluatableExpr* evaluatableArg() const 
          {return evaluatableChild(0);}

          /** */
          virtual Set<MultiIndex> argMultiIndices(const Set<MultiIndex>& multiIndices) const 
          {return multiIndices;}
          

          /** */
          virtual Set<MultiSet<int> > 
          argActiveFuncs(const Set<MultiSet<int> >& activeFuncID) const 
          {return activeFuncID;}
          
          /** */
          void addActiveFuncs(const EvalContext& context,
                              const Set<MultiSet<int> >& activeFuncIDs) const ;
          
          /** */
          const Set<MultiSet<int> >& getActiveFuncs(const EvalContext& context) const ;
        private:
          mutable Map<EvalContext, Set<MultiSet<int> > > allActiveFuncs_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
