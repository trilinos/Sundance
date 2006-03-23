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

#ifndef SUNDANCE_NONLINEARUNARYOP_H
#define SUNDANCE_NONLINEARUNARYOP_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceNonlinearUnaryOpEvaluator.hpp"
#include "SundanceNonlinearExpr.hpp"



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
    class NonlinearUnaryOp : public UnaryExpr,
                             public NonlinearExpr,
                             public GenericEvaluatorFactory<NonlinearUnaryOp, NonlinearUnaryOpEvaluator>
    {
    public:
      /** construct with an argument and the functor defining the operation */
      NonlinearUnaryOp(const RefCountPtr<ScalarExpr>& arg, 
                       const RefCountPtr<UnaryFunctor>& op);

      /** virtual destructor */
      virtual ~NonlinearUnaryOp() {;}

      /** Preprocessing step to determine which functional 
       * derivatives are nonzero */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const ;

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}


      /** Access to the operator */
      const UnaryFunctor* op() const {return op_.get();}

      /** 
       * Given a set of active function combinations, get the active
       * function combinations requested of the argument
       */
      virtual Set<MultiSet<int> > 
      argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs,
                     int maxOrder) const ;

      /** Given a set of required multiindices, find the set of multiindices
       * required of the argument */
      Set<MultiIndex> argMultiIndices(const Set<MultiIndex>& miSet) const ;
    private:

      RefCountPtr<UnaryFunctor> op_;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
