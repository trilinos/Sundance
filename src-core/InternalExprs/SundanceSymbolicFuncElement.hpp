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

#ifndef SUNDANCE_SYMBOLICFUNCELEMENT_H
#define SUNDANCE_SYMBOLICFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    class DiscreteFuncElement;
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * SymbolicFuncElement represents a scalar-valued element of a (possibly)
     * list-valued SymbolicFunction. 
     */
    class SymbolicFuncElement : public FuncElementBase,
                                virtual public EvaluatableExpr,
                                public GenericEvaluatorFactory<SymbolicFuncElement, SymbolicFuncElementEvaluator>
    {
    public:
      /** */
      SymbolicFuncElement(const string& name, const string& suffix,
                          int myIndex);
      
      /** virtual destructor */
      virtual ~SymbolicFuncElement() {;}

      /** Get my index into the master's list of elements */
      int myIndex() const {return myIndex_;}


      /** Specify that expressions involving this function are to be evaluated
       * with this function set to zero. Test functions should always be
       * evaluated at zero. For unknown functions, 
       * substituting zero is appropriate for computing
       * the functional derivatives that arise in a linear problem.
       * */
      void substituteZero() const ;

      /** Specify that expressions involving this function are to be evaluated
       * with this function set to the discrete function (or constant parameter) \f$u_0\f$. 
       * This is appropriate for computing
       * the functional derivatives that arise in a nonlinear expression
       * being linearized about \f$u_0\f$. 
       */
      void substituteFunction(const RefCountPtr<DiscreteFuncElement>& u0) const ;

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      const EvaluatableExpr* evalPt() const {return evalPt_.get();}

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      EvaluatableExpr* evalPt() {return evalPt_.get();}


      /** */
      bool evalPtIsZero() const ;


      /** \name Preprocessing */
      //@{
      /** 
       * Determine which functional and spatial derivatives are nonzero in the
       * given context. We also keep track of which functional derivatives
       * are known to be constant, which can simplify evaluation. 
       */
      virtual void findNonzeros(const EvalContext& context,
                                const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                bool regardFuncsAsConstant) const ;

          
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
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      mutable RefCountPtr<EvaluatableExpr> evalPt_;

      mutable Array<int> evalPtDerivSetIndices_;

      int myIndex_;
      

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
