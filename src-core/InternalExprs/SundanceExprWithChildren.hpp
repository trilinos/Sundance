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

          /** Append to the set of func IDs present in 
           * this expression. */
          virtual void accumulateFuncSet(Set<int>& funcIDs,
                                         const Set<int>& activeFuncs) const ;

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
                                    const Set<MultiSet<int> >& activeFuncIDs,
                                    bool regardFuncsAsConstant) const ;
          
          /** Return true if any child returns true. The sum expression
           * will override this requiring all children to return true */
          virtual bool allTermsHaveTestFunctions() const ;

          /** Test whether this expr contains a test function. 
           * If any child contains a test, return true. */
          virtual bool hasTestFunctions() const ;

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
