/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALUATABLEEXPR_H
#define SUNDANCE_EVALUATABLEEXPR_H



#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceMap.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalVectorArray.hpp"
#include "SundanceSparsityPattern.hpp"
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
     * Evaluatable expression subtypes may use multiple inheritance.
     * EvaluatableExpr contains many data structures that should
     * only exist once per object, so EvaluatableExpr is a virtual
     * base class. Classes deriving from EvaluatableExpr should use
     * virtual public inheritance.
     *
     * Evaluation is a logically const operation, and preprocessing
     * is an implementation detail, so most evaluation and preprocessing
     * related methods of EvaluatableExpr are declared const, and
     * the internal data is declared mutable.
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
     * <h3> Description of sparsity structure </h3>
     *
     * Information about exactly which derivatives are nonzero
     * and required
     * at a given node in a given region is
     * stored in a SparsityPattern object. Since the derivatives required
     * will depend on region, there will be a SparsityPattern
     * object for each region, or more precisely for each set of multiple
     * derivatives required at the root.
     *
     * A SparsityPattern object describes the set of all derivatives
     * that will be nonzero during the course of a particular
     * root-level evaluate() call. At nodes farther down the graph,
     * some of those derivatives may become zero. To simplify the
     * bookkeeping by preserving consistent numbering between many
     * nodes, the SparsityPattern object still includes those
     * derivatives but marks them with a ZeroDeriv value of a
     * DerivState enum. Other derivatives may be identifiable as
     * constant at a node, in which case they are marked with a
     * ConstantDeriv state. Derivatives that are nonzero and
     * non-constant are marked with a VectorDeriv state.
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
     *
     *
     * <h4> Managing multiple calls during evaluation </h4>
     *
     * Since expressions have shallow copy behavior,
     * any expression object may
     * appear multiple times in an expression graph, and thus during an
     * evaluation of that graph some node may have their evaluate()
     * methods called multiple times. Clearly it is wasteful to
     * repeatedly evaluate the same expression, so a system for
     * avoiding multiple evaluations is required.
     *
     * Consider the sequence of expressions \code Expr d = a + b;
     * Expr e = a + c*d + sqrt(d); \endcode Upon evaluation of
     * <tt>e</tt>, there will clearly be one call to c.evaluate()
     * and two to d.evaluate(). But how many calls are there to
     * a.evaluate() and b.evaluate()? Both are called once in the
     * evaluation of d, but then a.evaluate() is called one more
     * time in the evaluation of e. However, since we will cache d's
     * value, any subsequent calls to d.evaluate() will not
     * recursively call a.evaluate(). So the number of evaluate()
     * calls for each of \f$\{a,b,c,d\}\f$ is \f$\{2,1,1,2\}\f$.
     *
     * The number of calls to each node's evaluate during an
     * evaluation of the root node is the number of edges
     * entering the node from above in the expression graph.
     * This can be computed as in the following pseudocode:
     * \code
     * 1. Initialize all nodes' reference counts to zero
     * 2. Call root->resetReferenceCount(), where
     *
     * void resetReferenceCount()
     * {
     *   increment refCount;
     *   if (refCount == 1)
     *   {
     *     for all children {child->resetReferenceCount();}
     *   }
     * }
     * \endcode
     * During a subsequent evaluation, if a node's reference count
     * is greater than zeo then its value should be copied into a cache
     * for further use. Each call to evaluate() decrements the reference
     * count, and when the count reaches zero the storage allocated
     * for the cached values can be placed back on a stack and
     * later reused.
     */
    class EvaluatableExpr : public virtual ScalarExpr,
                            public TSFExtended::ObjectWithVerbosity<EvaluatableExpr>
    {
    public:
      /** */
      EvaluatableExpr();

      /** */
      virtual ~EvaluatableExpr(){;}

      /**
       * Indicate whether the given functional derivative is nonzero.
       * Each subclass implementing the EvaluatableExpr interface
       * must use the appropriate rules of calculus to make
       * that determination for the operation or object it represents.
       */
      virtual bool hasNonzeroDeriv(const MultipleDeriv& f) const=0;

      /**
       * Find all functions and their derivatives appearing
       * beneath my level in the tree, regardless of whether the
       * derivative wrt that function will prove to be
       * structurally zero. This method does not do any math to
       * check which lower-order derivatives propagate up the
       * tree, and so might generate false positives. This is just
       * a rough means of finding the set of functions and
       * derivatives that need to be checked in the chain rule
       * application in DiffOp::hasNonzeroDeriv().
       */
      virtual void getRoughDependencies(Set<Deriv>& funcs) const = 0 ;
      
      /** Return the set of this expression's 
       * nonzero functional derivatives through second order */
      DerivSet identifyNonzeroDerivs() const ;

      /**
       * Do preprocessing for efficient evaluation in the given region, and 
       * then create and install the appropriate evaluator as prescribed
       * by the evaluator factory.
       */
      virtual int setupEval(const RegionQuadCombo& region,
                            const EvaluatorFactory* factory) const = 0;

      /**
       * Evaluate this expression in the given region, putting the results
       * of the evaluation in the results argument. 
       */
      void evaluate(const EvalManager& mgr,
                    RefCountPtr<EvalVectorArray>& results) const ;


      /** Look up the index at which the sparsity information for this
       * region is stored. */
      int getDerivSetIndex(const RegionQuadCombo& region) const ;


      /** Return the sparsity pattern for the given deriv set */
      const RefCountPtr<SparsityPattern>& sparsity(int derivSetIndex) const
      {return sparsityPatterns_[derivSetIndex];}

      /** */
      bool hasWorkspace() const ;


      /** \name Reference counting */
      //@{
      /** Return the reference count for evaluations. This is the
       * number of times this expression will be evaluated
       * before its storage can be released. */
      int& refCount() const {return evalRefCount_;}

      /**
       * Recompute the evaluation reference counts for this node and its
       * children.
       */
      virtual void resetReferenceCount() const {refCount()++;}
      //@}

      /** flush the cache of computed values */
      virtual void flushResultCache() const ;


      /** See if this region is new */
      bool checkForKnownRegion(const RegionQuadCombo& region,
                               bool& derivSetIsKnown) const ;

      /** Create entries in the results array for the all the
       * elements of this set of derivatives, and insert the
       * list of indices into derivSetIndexToResultIndicesMap_.
       * Concurrently create a SparsityPattern object and an Evaluator.
       * 
       * @return the index of the newly-created list of indices in
       * derivSetIndexToResultIndicesMap_
       */
      int registerRegion(const RegionQuadCombo& region,
                         bool derivSetIsKnown,
                         const DerivSet& derivs,
                         const EvaluatorFactory* factory) const ;

      /** Indicate whether this region has been set up */
      bool knowsRegion(const RegionQuadCombo& region) const
      {return regionToDerivSetIndexMap_.containsKey(region);}


      /** Return the deriv set stored at the given index */
      const DerivSet& getDerivSet(int derivSetIndex) const 
      {return derivSets_[derivSetIndex];}

      /** \name Determination of Deriv superset */
      //@{
      /** */
      virtual void resetDerivSuperset() const {currentDerivSuperset() = DerivSet();}
      /** */
      virtual void findDerivSuperset(const DerivSet& derivs) const ;
      //@}

      /** */
      static Time& evalTimer() 
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Expr evaluation"); 
        return *rtn;
      }
      /** */
      static Time& derivIDTimer() 
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Expr derivID"); 
        return *rtn;
      }

      /** */
      static const EvaluatableExpr* getEvalExpr(const Expr& expr);

    protected:

      /** \name Caching of hasNonzeroDeriv() results */
      //@{
      /** Indicate whether the cache contains a result for the given
       * derivative. */
      bool derivHasBeenCached(const MultipleDeriv& d) const
      {return derivNonzeronessCache_.containsKey(d);}

      /** Look up the nonzeroness of the given derivative in the cache. */
      bool getCachedDerivNonzeroness(const MultipleDeriv& d) const
      {return derivNonzeronessCache_.get(d);}

      /**
       * Add to the cache the nonzeroness of the given derivative
       */
      void addDerivToCache(const MultipleDeriv& d, bool hasF) const
      {derivNonzeronessCache_.put(d, hasF);}
      //@}

      /** Return the evaluator to be used for the given deriv set */
      const Evaluator* evaluator(int derivSetIndex) const 
      {return evaluators_[derivSetIndex].get();}

      /** */
      DerivSet& currentDerivSuperset() const {return currentDerivSuperset_;}

      /** */
      bool& resultCacheIsValid() const {return resultCacheIsValid_;}

    private:
      /**
       * Map from evaluation region to index in array of derivative
       * sets.
       *
       * The mapping from region to derivative set is possibly
       * many-to-one, because two or more regions may require
       * the same set of derivatives from this node. Rather than replicate
       * identical derivative sets in
       * a <tt> Map<RegionQuadCombo, DerivSet ></tt> object,
       * create an array of derivative sets and map regions to
       * that array index.
       */
      mutable Map<RegionQuadCombo, int> regionToDerivSetIndexMap_;

      /**
       * Derivative sets are usually referred to by index. This map
       * lets us obtain the index in the deriv set list
       * at which a given deriv set is stored.
       */
      mutable Map<DerivSet, int> derivSetToDerivSetIndexMap_;

      /** 
       *
       */
      mutable Array<DerivSet> derivSets_;

      /**
       * Array of sparsity patterns, indexed by deriv set.
       */
      mutable Array<RefCountPtr<SparsityPattern> > sparsityPatterns_;

      /** 
       * Array of evaluators, indexed by deriv set.
       */
      mutable Array<RefCountPtr<Evaluator> > evaluators_;

      /**
       * Determination of whether a given derivative is nonzero
       * is expensive and need only be done once per node. This data
       * member stores for each derivative that has been processed
       * a bool indicating whether its derivative is nonzero.
       */
      mutable Map<MultipleDeriv, bool> derivNonzeronessCache_;

      /**
       * The number of times this node needs to be evaluated in the
       * current state of the evaluation.
       */
      mutable int evalRefCount_;

      /** 
       *
       */
      mutable bool resultCacheIsValid_;

      /**
       * Cache of results so that previously computed values can be
       * reused. 
       */
      mutable RefCountPtr<EvalVectorArray> resultCache_;

      /** */
      mutable DerivSet currentDerivSuperset_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
