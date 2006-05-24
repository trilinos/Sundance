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

#ifndef SUNDANCE_DIFFOP_H
#define SUNDANCE_DIFFOP_H

#include "SundanceDefs.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceDiffOpEvaluator.hpp"


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
    class DiffOp : public UnaryExpr
    {
    public:
      /** ctor */
      DiffOp(const MultiIndex& op, const RefCountPtr<ScalarExpr>& arg);

      /** virtual destructor */
      virtual ~DiffOp() {;}

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;


      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const ;


      /** */
      virtual Set<MultipleDeriv> 
      internalFindV(int order, const EvalContext& context) const ;


      /** */
      virtual Set<MultipleDeriv> 
      internalFindC(int order, const EvalContext& context) const ;


      /** */
      virtual RefCountPtr<Array<Set<MultipleDeriv> > > 
      internalDetermineR(const EvalContext& context,
                         const Array<Set<MultipleDeriv> >& RInput) const ;
      
      
      
      /** */
      const Deriv& myCoordDeriv() const {return myCoordDeriv_;}

      /** */
      const MultiIndex& mi() const {return mi_;}

      /** Get the functions that are required in the evaluation
       * of the multiple deriv d */
      const SundanceUtils::Set<Deriv>& requiredFunctions(const MultipleDeriv& d) const 
      {return requiredFunctions_[d];}

      /** */
      bool requiresFunctionsToEval(const MultipleDeriv& d) const 
      {return requiredFunctions_.containsKey(d);}

    
      
      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}


      /** Determine the functional derivatives of this expression that
       * are produced through the chain rule
       * by the presence of the given derivative of the argument. 
       */
      void getResultDerivs(const MultipleDeriv& argDeriv,
                           const DerivState& sourceState,
                           Map<MultipleDeriv, DerivState>& isolatedTerms,
                           Map<MultipleDeriv, Deriv>& funcTerms) const ;


      /** */
      virtual Set<MultiIndex> argMultiIndices(const Set<MultiIndex>& myMultiindices) const ;
        
      /** */
      bool ignoreFuncTerms() const {return ignoreFuncTerms_;}

      /** */
      virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                         const EvalContext& context) const ;

      /** 
       * Given a set of active function combinations, get the active
       * funcs for the next higher order of differentiation.
       */
      virtual Set<MultiSet<int> > 
      argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs,
                     int maxOrder) const ;



    private:

      
      bool canBackOutDeriv(const Deriv& d, const MultiIndex& x, 
                           Deriv& rtnDeriv) const ;

      

      MultiIndex mi_;


      Deriv myCoordDeriv_;

      mutable Map<MultipleDeriv, SundanceUtils::Set<Deriv>, 
                  increasingOrder<MultipleDeriv> > requiredFunctions_;


      mutable bool ignoreFuncTerms_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
