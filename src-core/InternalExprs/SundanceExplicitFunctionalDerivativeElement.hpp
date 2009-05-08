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

#ifndef SUNDANCE_EXPLICITFUNCTIONALDERIVATIVEELEMENT_H
#define SUNDANCE_EXPLICITFUNCTIONALDERIVATIVEELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceExprWithChildren.hpp"


namespace SundanceCore
{
using namespace SundanceUtils;
using namespace Teuchos;

using std::string;
using std::ostream;


/** 
 * Scalar element of an explicit user-level functional derivative
 */
class ExplicitFunctionalDerivativeElement 
  : virtual public UnaryExpr
{
public:
  /** */
  ExplicitFunctionalDerivativeElement(
    const RefCountPtr<ScalarExpr>& arg,
    const Deriv& fd
    );
            
  /** virtual destructor */
  virtual ~ExplicitFunctionalDerivativeElement() {;}

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


  /** Write self in text form */
  virtual ostream& toText(ostream& os, bool paren) const ;

  /** Write self in Latex form */
  virtual ostream& toLatex(ostream& os, bool paren) const ;

  /** Write in XML */
  virtual XMLObject toXML() const ;

  /** */
  virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}

  /** */
  void reset() const ;

  /** */
  Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

  /** */
  virtual bool lessThan(const ScalarExpr* other) const ;


  /** */
  const Deriv& fd() const {return fd_;}

protected:
private:
  Deriv fd_;
};

}


#endif
