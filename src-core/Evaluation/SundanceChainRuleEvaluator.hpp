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

#ifndef SUNDANCE_CHAINRULEEVALUATOR_H
#define SUNDANCE_CHAINRULEEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceExprWithChildren.hpp"
#include "SundanceChainRuleSum.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore 
{
  namespace Internal 
  {
    /** 
     *
     */
    class ChainRuleEvaluator : public SubtypeEvaluator<ExprWithChildren>
    {
    public:

      /** */
      ChainRuleEvaluator(const ExprWithChildren* expr, 
                         const EvalContext& context);

      /** */
      virtual ~ChainRuleEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                                Array<double>& constantResults,
                                Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      int numChildren() const {return childEvaluators_.size();}

      /** */
      int constArgDerivIndex(const MultiSet<int>& df) const ;

      /** */
      int varArgDerivIndex(const MultiSet<int>& df) const ;

      /** */
      TEUCHOS_TIMER(chainRuleEvalTimer, "chain rule evaluation");


      /** */
      const Array<Array<int> >& nComps(int N, int n) const ;

      /** */
      void resetNumCalls() const ;

      /** */
      virtual void evalArgDerivs(const EvalManager& mgr,
                                 const Array<RefCountPtr<Array<double> > >& constArgRes,
                                 const Array<RefCountPtr<Array<RefCountPtr<EvalVector> > > >& vArgResults,
                                 Array<double>& constArgDerivs,
                                 Array<RefCountPtr<EvalVector> >& varArgDerivs) const = 0 ;


      static Set<MultiSet<MultipleDeriv> > chainRuleBins(const MultipleDeriv& d,
                                                      const MultiSet<int>& q);
      
    protected:
      /** The init() function should be called from the derived class ctors */
      void init(const ExprWithChildren* expr, 
                const EvalContext& context);

      /** */
      void addConstArgDeriv(const MultiSet<int>& df, int index);

      /** */
      void addVarArgDeriv(const MultiSet<int>& df, int index);

      /** */
      const Evaluator* childEvaluator(int i) const {return childEvaluators_[i].get();}

      /** */
      const SparsitySuperset* childSparsity(int i) const {return childSparsity_[i].get();}

      static MultipleDeriv makeMD(const Array<Deriv>& d) ;

      /** Returns the binomial coefficient */
      double choose(int N, int n) const ;
      
      /** Returns the factorial of n */
      double fact(int n) const ;


      /** Returns the stirling number of the second kind */
      double stirling2(int n, int k) const ;

      /** */
      int derivComboMultiplicity(const MultiSet<MultipleDeriv>& b) const ;


    private:

      
      Array<RefCountPtr<ChainRuleSum> > expansions_;

      Array<RefCountPtr<Evaluator> > childEvaluators_;

      Array<RefCountPtr<SparsitySuperset> > childSparsity_;

      Map<MultiSet<int>, int> constArgDerivMap_;

      Map<MultiSet<int>, int> varArgDerivMap_;

      int zerothDerivResultIndex_;

      bool zerothDerivIsConstant_;

      static Map<OrderedPair<int, int>, Array<Array<int> > >& compMap() ;
    };

    /** */
    MultipleDeriv makeDeriv(const Expr& a);

    /** */
    MultipleDeriv makeDeriv(const Expr& a, const Expr& b);

    /** */
    MultipleDeriv makeDeriv(const Expr& a, const Expr& b, const Expr& c);

    
  }
}
               
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif