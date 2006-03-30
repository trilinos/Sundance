/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
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
     * determination of the sparsity pattern of the matrix of functional
     * derivatives, and managing Evaluator objects. Because we will 
     * need to evaluate expressions in different ways depending
     * on context, we shouldn't do evaluation in a method
     * of EvaluatableExpr, but rather delegate evaluation to a
     * separate class, Evaluator. EvaluatableExpr stores a suite of
     * Evaluator objects, indexed by EvalCOntext.

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
                                                 const Set<MultiSet<int> >& activeFuncIDs,
                                                 bool failIfNotFound) const ;



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
      static unsigned int maxFuncDiffOrder() {static int rtn=3; return rtn;}

      /** */
      static RefCountPtr<Set<int> > getFuncIDSet(const Expr& funcs);

      /** Filter the input active function list to eliminate functions
       * not appearing in the list of functional dependencies */
      Set<MultiSet<int> > filterActiveFuncs(const Set<MultiSet<int> >& inputActiveFuncs) const ;


      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const = 0 ;

      /** */
      const Set<MultipleDeriv>& findW(int order, 
                                      const EvalContext& context) const ;
      /** */
      Set<MultipleDeriv> setProduct(const Set<MultipleDeriv>& a,
                                    const Set<MultipleDeriv>& b) const ;
      
      /** */
      void determineR(const EvalContext& context,
                      const Array<Set<MultipleDeriv> >& RInput) const ;

      /** */
      virtual RefCountPtr<Array<Set<MultipleDeriv> > > 
      internalDetermineR(const EvalContext& context,
                         const Array<Set<MultipleDeriv> >& RInput) const ;

      /** */
      const Set<MultipleDeriv>& getR(int order, const EvalContext& context) const ;

      
      /** */
      Array<Set<MultipleDeriv> > 
      computeInputR(const EvalContext& context,
                    const Array<Set<MultiSet<int> > >& funcIDCombinations,
                    const Array<Set<MultiIndex> >& spatialDerivs) const ;
      
                
      
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



      /** 
       * Computes
       * \f[ 
       * T_x(S) = \left\{\mu \vert \lambda\in S \wedge (\alpha_\lambda + x=\alpha_\mu)
       *           \wedge (u_\lambda = u_\mu) \right\}
       * \f]
       */
      static Set<MultipleDeriv> computeT(int xDir, const Set<MultipleDeriv>& S);

      /** 
       * Computes
       * \f[ 
       * K_x(S) = \left\{\mu \vert \mu\in S \wedge p_W(D_{\alpha_\mu + x} u_\mu) \right\} 
       * \f]
       */
      static Set<MultipleDeriv> computeK(int xDir, const Set<MultipleDeriv>& S);

      /** */
      static Set<MultipleDeriv> computeH(bool pred, const Set<MultipleDeriv>& S);
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

      mutable Array<Map<EvalContext, Set<MultipleDeriv> > > contextToRTableMap_;
      mutable Array<Map<EvalContext, Set<MultipleDeriv> > > contextToWTableMap_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
