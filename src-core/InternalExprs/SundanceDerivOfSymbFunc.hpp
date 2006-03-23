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

#ifndef SUNDANCE_DERIVOFSYMBFUNC_H
#define SUNDANCE_DERIVOFSYMBFUNC_H

#include "SundanceDefs.hpp"
#include "SundanceDiffOp.hpp"
#include "SundanceDerivOfSymbFuncEvaluator.hpp"


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
     * Specialization of DiffOp to the case where the argument is a 
     * symbolic function, allowing optimized evaluation.
     */
    class DerivOfSymbFunc : public DiffOp
    {
    public:
      /** ctor */
      DerivOfSymbFunc(const MultiIndex& op, 
                      const RefCountPtr<ScalarExpr>& arg);

      /** virtual destructor */
      virtual ~DerivOfSymbFunc() {;}

      /** */
      int funcID() const {return funcID_;}

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

      /** */
      virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                         const EvalContext& context) const ;

      /** */
      virtual Set<MultiSet<int> > argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs,
                                                 int maxOrder) const ;

      /** */
      FunctionalDeriv* representMeAsFunctionalDeriv() const ;

    private:
      int funcID_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
