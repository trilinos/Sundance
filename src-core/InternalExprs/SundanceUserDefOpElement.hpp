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

#ifndef SUNDANCE_USERDEFOPELEMENT_H
#define SUNDANCE_USERDEFOPELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefFunctorElement.hpp"
#include "SundanceUserDefOpEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;


  class UserDefOp;

  namespace Internal
    {
      /** 
       * Scalar element of a vector-valued user-defined expression.
       */
      class UserDefOpElement : virtual public ExprWithChildren
        {
        public:
          /** */
          UserDefOpElement(const Array<RefCountPtr<ScalarExpr> >& args,
            const RefCountPtr<SundanceUtils::Map<EvalContext, RefCountPtr<const UserDefOpCommonEvaluator> > >& evalMap,
            const RefCountPtr<const UserDefFunctorElement>& functorElement);

          /** virtual destructor */
          virtual ~UserDefOpElement() {;}

          /** Return the index of this element into 
           * the list-valued user defined op */
          int myIndex() const {return functorElement_->myIndex();}


          /** Write self in text form */
          virtual ostream& toText(ostream& os, bool paren) const ;

          /** Write self in Latex form */
          virtual ostream& toLatex(ostream& os, bool paren) const ;

          /** Write in XML */
          virtual XMLObject toXML() const ;

          /** */
          virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

          /** Access to the functor underlying this object */
          const UserDefFunctorElement* functorElement() const 
          {return functorElement_.get();}

          /** */
          void reset() const ;

          /** */
          Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                     const EvalContext& context) const ;

          /** */
          virtual void getArgDerivIndices(const Array<int>& orders,
                                          SundanceUtils::Map<MultiSet<int>, int>& varArgDerivs,
                                          SundanceUtils::Map<MultiSet<int>, int>& constArgDerivs) const ;

        protected:
          /** Get an evaluator that will be common to all vector elements of 
           * this operator */
          RefCountPtr<const UserDefOpCommonEvaluator> 
          getCommonEvaluator(const EvalContext& context) const ;
          
        private:
          mutable RefCountPtr<SundanceUtils::Map<EvalContext, RefCountPtr<const UserDefOpCommonEvaluator> > > commonEvaluatorsMap_;
          const RefCountPtr<const UserDefFunctorElement> functorElement_;
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif