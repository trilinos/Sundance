/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXPRWITHCHILDREN_H
#define SUNDANCE_EXPRWITHCHILDREN_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceScalarExpr.hpp"
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
       * ExprWithChildren is a base class for any evaluatable expression
       * that has child nodes, for example, sums and unary operators.
       * ExprWithChildren adds nothing new to the expr interface, but 
       * provides some common utilities for getting children
       * and recursing to children.
       */
      class ExprWithChildren : public virtual EvaluatableExpr
        {
        public:
          /** construct with a list of child operands */
          ExprWithChildren(const Array<RefCountPtr<ScalarExpr> >& children);

          /** virtual dtor */
          virtual ~ExprWithChildren() {;}

          /**
           * Do preprocessing to set up sparse evaluation in the given region 
           */
          virtual void setupEval(const EvalContext& context) const ;

          /** Determine whether this expression is constant. It will
           * be constant if all children are constant. */
          virtual bool isConstant() const ;

          /** Append to the set of unknown func IDs present in 
           * this expression. */
          virtual void accumulateUnkSet(Set<int>& unkIDs) const ;

          /** Append to the set of test func IDs present in this 
           * expression. */
          virtual void accumulateTestSet(Set<int>& testIDs) const ;

          /** Return the number of children */
          int numChildren() const {return children_.size();}
          
          /** downcast the i-th to an evaluatable expr */
          const EvaluatableExpr* evaluatableChild(int i) const ;

          /** downcast the i-th to a scalar expr */
          const ScalarExpr* scalarChild(int i) const 
          {return children_[i].get();}

          /** Get a handle to the i-th child */
          Expr child(int i) const {return Expr::handle(children_[i]);}

          /** 
           * Generic ExprWithChildren objects find the sparsity to be
           * the union of the children's sparsity patterns. DiffOp and
           * ProductExpr will need to override this method.
           */
          virtual void findNonzeros(const EvalContext& context,
                                    const Set<MultiIndex>& multiIndices,
                                    bool regardFuncsAsConstant) const ;
          
          /** Return true if any child returns true. The sum expression
           * will override this requiring all children to return true */
          virtual bool allTermsHaveTestFunctions() const ;

          /** */
          virtual void showSparsity(ostream& os, 
                                    const EvalContext& context) const ;

          /** */
          virtual void getUnknowns(Set<int>& unkID, Array<Expr>& unks) const ;

          
          /** */
          virtual int countNodes() const ;
        private:
          Array<RefCountPtr<ScalarExpr> > children_;
      };          

    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
