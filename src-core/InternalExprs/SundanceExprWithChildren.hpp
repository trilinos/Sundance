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
      /** */
      class ExprWithChildren : public virtual EvaluatableExpr
        {
        public:
          /** construct with left and right operands */
          ExprWithChildren(const Array<RefCountPtr<ScalarExpr> >& children);

          /** virtual dtor */
          virtual ~ExprWithChildren() {;}

          /** */
          virtual bool isConstant() const ;

          /** Append to the set of unknown func IDs present in 
           * this expression. */
          virtual void accumulateUnkSet(Set<int>& unkIDs) const ;

          /** Append to the set of test func IDs present in this 
           * expression. */
          virtual void accumulateTestSet(Set<int>& testIDs) const ;

          
          /** downcast the i-th to an evaluatable expr */
          const EvaluatableExpr* evaluatableChild(int i) const ;

          /** downcast the i-th to a scalar expr */
          const ScalarExpr* scalarChild(int i) const 
          {return children_[i].get();}

          /** Get a handle to the i-th child */
          Expr child(int i) const {return Expr::handle(children_[i]);}


          /** Return the index by which the given deriv set is
           * known by the given child */
          int childDerivSetIndex(int child, int derivSetIndex) const 
          {return childDerivSetIndices_[child][derivSetIndex];}

          /**
           * Find all functions and their derivatives beneath my level
           * in the tree. 
           */
          virtual void getRoughDependencies(Set<Deriv>& funcs) const ;

          /**
           * Recompute the reference counts for this node and its
           * children.
           */
          virtual void resetReferenceCount() const ;

          /** flush the cache of computed values */
          virtual void flushResultCache() const ;

          /** */
          virtual bool hasWorkspace() const ;

          /** */
          virtual void findDerivSuperset(const DerivSet& derivs) const ;

          /** */
          virtual void resetDerivSuperset() const ;

          /** Do the steps that are common over all binary ops
           * for setting up evaluation of a new deriv set.
           * @return index of the deriv set
           */
          virtual int setupEval(const EvalContext& region,
                                const EvaluatorFactory* factory,
                                bool regardFuncsAsConstant) const ;

          /** Return the set of derivatives required by the operands
           * of this expression given that this expression
           * requires the set d. For all expressions other than DiffOp,
           * the operand derivative set is identical to the input derivative
           * set. DiffOp will require a different set of derivatives from
           * its operand. */
          virtual Array<DerivSet>
          derivsRequiredFromOperands(const DerivSet& d) const ;

        private:
          Array<RefCountPtr<ScalarExpr> > children_;
          
          mutable Array< Array<int> > childDerivSetIndices_;
        };
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
