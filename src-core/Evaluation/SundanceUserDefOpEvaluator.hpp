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


#ifndef SUNDANCE_USERDEFOPEVALUATOR_H
#define SUNDANCE_USERDEFOPEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore 
{
  class UserDefOp;

  namespace Internal 
  {

    class SymbolicFuncElementEvaluator;
    
    /**
     *
     */
    class UserDefOpEvaluator : public SubtypeEvaluator<UserDefOp>
    {
    public:
      /** */
      UserDefOpEvaluator(const UserDefOp* expr,
                         const EvalContext& context);

      /** */
      virtual ~UserDefOpEvaluator(){;}

      /** */
      virtual void internalEval(const EvalManager& mgr,
                                Array<double>& constantResults,
                                Array<RefCountPtr<EvalVector> >& vectorResults) const ;

      /** */
      virtual void resetNumCalls() const ;

      /** */
      TEUCHOS_TIMER(evalTimer, "user defined nonlinear op evaluation");

    protected:
      const RefCountPtr<SparsitySuperset>& childSparsity(int i) const
      {return childSparsity_[i];}

      const EvaluatableExpr* childExpr(int i) const
      {return childExpr_[i];}

      const RefCountPtr<Evaluator>& childEval(int i) const
      {return childEval_[i];}

      void evalChildren(const EvalManager& mgr,
                        Array<Array<double> >& constResults,
                        Array<Array<RefCountPtr<EvalVector> > >& vecResults) const ;

      void evalOperator(int numPoints,
                        const Array<double>& constantArg,
                        const Array<RefCountPtr<EvalVector> >& vectorArg,
                        const Array<int>& constantArgPtr,
                        const Array<int>& vectorArgPtr,
                        RefCountPtr<EvalVector>& opResults) const ;

      
    private:
      Array<const EvaluatableExpr*> childExpr_;

      Array<RefCountPtr<SparsitySuperset> > childSparsity_;

      Array<RefCountPtr<Evaluator> > childEval_;

      int maxOrder_;
      int d0ResultIndex_;
      Array<int> d0ArgDerivIndex_;
      Array<int> d0ArgDerivIsConstant_;
      Array<int> constantArgPtr_;
      Array<int> vectorArgPtr_;
    }; 
  }
}
                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  


#endif
