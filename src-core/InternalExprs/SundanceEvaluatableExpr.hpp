/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALUATABLEEXPR_H
#define SUNDANCE_EVALUATABLEEXPR_H



#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceMap.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore::Internal;

  using std::string;
  using std::ostream;
  class Expr;

  namespace Internal
  {
    class MultipleDeriv;
    class EvalManager;
    class Evaluator;
    class EvaluatorFactory;

    /**
     * Class EvaluatableExpr provides common functionality for
     * management of caching evaluation results, preprocessing
     * to determine the sparsity pattern of the matrix of functional
     * derivatives, and part of the interface for evaluation of expressions.
     *
     * Evaluation is a logically const operation, and preprocessing
     * is an implementation detail, so most evaluation and preprocessing
     * related methods of EvaluatableExpr are declared const, and
     * the internal data is declared mutable.
     *
     * <h3> Evaluation context </h3>
     *
     * A given expression might be evaluated in a number of different
     * contexts, in which different spatial or functional derivatives
     * are required. Naturally, we will never want to evaluate quantities
     * that are not required in the current evaluation context, so we will
     * build sparsity tables and evaluation instructions on a 
     * context-by-context basis. 
     *
     * <h3> Sparsity determination </h3>
     *
     * It is important to know the sparsity structure of the matrix
     * of partial functional derivatives. To determine whether an
     * EvaluatableExpr* <tt>e</tt> has a nonzero partial derivative
     * <tt>d,</tt> call <tt>e->hasNonzeroDeriv(d).</tt>
     *
     * Since sparsity determination does take some time in a complicated
     * expression, the nonzeroness of a given derivative is cached
     * once determined. The cache object is derivNonzeronessCache_,
     * which is declared mutable to maintain the logical constness
     * of hasNonzeroDeriv().
     *
     * Information about exactly which derivatives are nonzero
     * and required
     * at a given node in a given region is
     * stored in a SparsitySuperset object. Since the derivatives required
     * will depend on region, there will be a SparsitySuperset
     * object for each region. 
     *
     * <h3> Evaluation </h3>
     *
     * The interface for evaluation is the evaluate() method. 
     * Arguments are
     * a const reference to an EvalManager and
     * a smart pointer to a EvalVectorArray. The EvalManager
     * EvalManager is to supply any provide callbacks for evaluation of objects
     * such as discrete functions which require knowledge of the
     * mesh and field data structures to be carried out. The
     * EvalVectorArray object contains vectors for the results of
     * all functional derivatives evaluated at this node.
     */

    class EvaluatableExpr : public virtual ScalarExpr,
                            public virtual EvaluatorFactory,
                            public TSFExtended::ObjectWithVerbosity<EvaluatableExpr>
    {
      typedef OrderedQuartet<EvalContext, 
                             Set<MultiIndex>,
                             Set<MultiSet<int> >,
                             bool> NonzeroSpecifier ;
    public:
      /** Ctor is empty, but has some internal initialization to do
       * and so must be called by all subclasses */
      EvaluatableExpr();

      /** virtual dtor */
      virtual ~EvaluatableExpr(){;}


      /** \name Evaluation */
      //@{
      /**
       * Evaluate this expression in the given region, putting the results
       * of the evaluation in the results argument. 
       */
      void evaluate(const EvalManager& mgr,
                    Array<double>& constantResults,
                    Array<RefCountPtr<EvalVector> >& vectorResults) const ;
      //@}

      /** \name Preprocessing */
      //@{
      /**
       * 
       */
      virtual void setupEval(const EvalContext& context) const ;

      /** Return the subset of nonzero derivatives required to evaluate
       * the given set of differential operators in the given context. */
      RefCountPtr<SparsitySubset> sparsitySubset(const EvalContext& context,
                                                 const Set<MultiIndex>& multiIndices,
                                                 const Set<MultiSet<int> >& activeFuncIDs) const ;



      /** Return the set of all nonzero derivatives
       * required in the given context */
      RefCountPtr<SparsitySuperset> sparsitySuperset(const EvalContext& context) const ;

      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const = 0 ;

      /** */
      bool nonzerosAreKnown(const EvalContext& context,
                            const Set<MultiIndex>& multiIndices,
                            const Set<MultiSet<int> >& activeFuncIDs,
                            bool regardFuncsAsConstant) const ;
          

      /** 
       * Return the polynomial order of this expr's dependence on the 
       * given spatial coordinate function. Non-polynomial functions
       * (e.g., sqrt()) are assigned order=-1. This convention is 
       * non-standard mathematically but allows us to use integer 
       * order variables for all functions. 
       * Unknown, test, and discrete functions are assumed non-polynomial.
       * By default, any otherwise unspecified expression is assumed 
       * to be non-polynomial. 
       */
      int orderOfSpatialDependency(int spatialDir) const 
      {return orderOfDependency_[spatialDir];}

      /** 
       * Return all function combinations for which the mixed
       * functional derivatives can be nonzero.
       */
      const Set<MultiSet<int> >& funcIDSet() const {return funcIDSet_;}

      /**
       * Return the set of functions that appear in this expression.
       */
      const Set<int>& funcDependencies() const {return funcDependencies_;}

      //@}
      
      /** \name Error checking */
      //@{
      /** Test whether all terms have test functions. We'll use this
       * to check the validity of weak forms */
      virtual bool allTermsHaveTestFunctions() const {return false;}
      //@}


      /** Utility to downcast an expression to an evaluatable expr. Throws
       * an exception if the cast fails. */
      static const EvaluatableExpr* getEvalExpr(const Expr& expr);

      /** Return the evaluator to be used for the given context */
      const RefCountPtr<Evaluator>& evaluator(const EvalContext& context) const; 
      /** */
      virtual void showSparsity(ostream& os, 
                                const EvalContext& context) const ;

      /** */
      virtual void getUnknowns(Set<int>& unkID, Array<Expr>& unks) const {;}

      /** */
      virtual int countNodes() const ;

      /** */
      virtual bool nodesHaveBeenCounted() const {return nodesHaveBeenCounted_;}

      /** */
      static int maxFuncDiffOrder() {static int rtn=3; return rtn;}

      /** */
      static RefCountPtr<Set<int> > getFuncIDSet(const Expr& funcs);

    protected:

      /** Record the evaluator to be used for the given context */
      void registerEvaluator(const EvalContext& context,
                             const RefCountPtr<Evaluator>& evaluator) const 
      {return evaluators_.put(context, evaluator);}

      /** Set the order of dependency of this expression on the
       * given spatial coordinate function. This method exists as a utility
       * to be called at construction time, and should probably not be
       * called at other times. */
      void setOrderOfDependency(int spatialDir, int order)
      {orderOfDependency_[spatialDir] = order;}

      /** */
      static bool isEvaluatable(const ExprBase* expr);

      /** */
      Map<EvalContext, RefCountPtr<Evaluator> >& evaluators() const 
      {return evaluators_;}


      /** */
      int maxOrder(const Set<MultiIndex>& m) const ;

      /** */
      void addKnownNonzero(const EvalContext& context,
                           const Set<MultiIndex>& multiIndices,
                           const Set<MultiSet<int> >& activeFuncIDs,
                           bool regardFuncsAsConstant) const ;


      /** */
      void addFuncIDCombo(const MultiSet<int>& funcIDSet);

      /** */
      void setFuncIDSet(const Set<MultiSet<int> >& funcIDSet);


    private:

      /** 
       * evaluators, indexed by context 
       */
      mutable Map<EvalContext, RefCountPtr<Evaluator> > evaluators_;

      /** 
       * supersets of nonzero derivatives to be computed, index by
       * context
       */
      mutable Map<EvalContext, RefCountPtr<SparsitySuperset> > sparsity_;

      /** Polynomial order of the dependency upon each coordinate direction */
      Array<int> orderOfDependency_;

      /** Set of function combinations appearing in nonzero mixed partials */ 
      Set<MultiSet<int> > funcIDSet_;

      /** */
      Set<int> funcDependencies_;

      mutable Set<NonzeroSpecifier> knownNonzeros_;

      mutable bool nodesHaveBeenCounted_; 
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
