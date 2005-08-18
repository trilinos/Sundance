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

#include "SundanceFunctionalPolynomial.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDerivOfSymbFunc.hpp"
#include "SundanceUnaryExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


FunctionalPolynomial::FunctionalPolynomial(const RefCountPtr<ScalarExpr>& expr)
  : funcs_(),
    coeffs_()
{
  TEST_FOR_EXCEPTION(!isConvertibleToPoly(expr.get()), InternalError, 
                     "FunctionalPolynomial ctor called with an argument that "
                     "cannot be converted to a polynomial");

  int funcID;
  MultiIndex mi;

  const SymbolicFuncElement* s 
    = dynamic_cast<const SymbolicFuncElement*>(expr.get());
  if (s != 0)
    {
      funcID = s->funcID();
      mi = MultiIndex();
    }
  
  const DerivOfSymbFunc* d
    = dynamic_cast<const DerivOfSymbFunc*>(expr.get());
  if (d != 0)
    {
      funcID = d->funcID();
      mi = d->mi();
    }

  if (d != 0 || s != 0)
    {
      MultiSet<OrderedPair<int, MultiIndex> > key;
      key.put(OrderedPair<int, MultiIndex>(funcID, mi));
      funcs_.put(funcID, expr);
      Expr coeff = 1.0;
      RefCountPtr<ScalarExpr> cPtr = rcp_dynamic_cast<ScalarExpr>(coeff.ptr());
      coeffs_.put(key, cPtr);
    }
}


FunctionalPolynomial::FunctionalPolynomial(const Map<int, RefCountPtr<ScalarExpr> >& funcs,
                                           const Map<MultiSet<OrderedPair<int, MultiIndex> >, 
                                           RefCountPtr<ScalarExpr> >& coeffs)
  : funcs_(funcs),
    coeffs_(coeffs)
{}

RefCountPtr<ScalarExpr> FunctionalPolynomial::
addPoly(const FunctionalPolynomial* other, int sign) const 
{
  Map<int, RefCountPtr<ScalarExpr> > funcs = funcs_;
  Map<MultiSet<OrderedPair<int, MultiIndex> >, RefCountPtr<ScalarExpr> > coeffs
    = coeffs_;

  
  for (Map<MultiSet<OrderedPair<int, MultiIndex> >, 
         RefCountPtr<ScalarExpr> >::const_iterator 
         i = other->coeffs_.begin(); i != other->coeffs_.end(); i++)
    {
      const MultiSet<OrderedPair<int, MultiIndex> >& key = i->first;
      Expr right = Expr::handle(i->second);
      Expr sum;
      
      if (coeffs.containsKey(key))
        {
          Expr left = Expr::handle(coeffs.get(key));

          if (sign > 0) sum = left + right;
          else sum = left - right;
        }
      else
        {
          if (sign > 0) sum = right;
          else sum = -right;
        }
      
      const ScalarExpr* se = dynamic_cast<const ScalarExpr*>(sum.ptr().get());
      TEST_FOR_EXCEPTION(se==0, InternalError,
                         "Sum could not be cast to scalar expression");
      coeffs.put(key, rcp_dynamic_cast<ScalarExpr>(sum.ptr()));
    }
  
  for (Map<int, RefCountPtr<ScalarExpr> >::const_iterator 
         i = other->funcs_.begin(); i != other->funcs_.end(); i++)
    {
      funcs.put(i->first, i->second);
    }

  Expr rtn = new FunctionalPolynomial(funcs, coeffs);
  return rcp_dynamic_cast<ScalarExpr>(rtn.ptr());
}


Evaluator* FunctionalPolynomial::createEvaluator(const EvaluatableExpr* expr,
                                                 const EvalContext& context) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError, "poly eval not ready");
  return (Evaluator*) 0;
}

void FunctionalPolynomial::findNonzeros(const EvalContext& context,
                                        const Set<MultiIndex>& multiIndices,
                                        const Set<MultiSet<int> >& activeFuncIDs,
                                        bool regardFuncsAsConstant) const
{
  TEST_FOR_EXCEPTION(true, RuntimeError, "poly eval not ready");
}


bool FunctionalPolynomial::isConvertibleToPoly(const ScalarExpr* expr)
{
  const SymbolicFuncElement* s 
    = dynamic_cast<const SymbolicFuncElement*>(expr);
  
  const DerivOfSymbFunc* d
    = dynamic_cast<const DerivOfSymbFunc*>(expr);

  const FunctionalPolynomial* p
    = dynamic_cast<const FunctionalPolynomial*>(expr);

  return (s != 0 || d != 0 || p != 0);
}


RefCountPtr<FunctionalPolynomial> FunctionalPolynomial::toPoly(const RefCountPtr<ScalarExpr>& expr)
{
  const FunctionalPolynomial* p
    = dynamic_cast<const FunctionalPolynomial*>(expr.get());

  if (p != 0) return rcp_dynamic_cast<FunctionalPolynomial>(expr);

  Expr rtn = new FunctionalPolynomial(expr);
  return rcp_dynamic_cast<FunctionalPolynomial>(rtn.ptr());
}


ostream& FunctionalPolynomial::toText(ostream& os, bool paren) const
{
  os << toXML();
  return os;
}

ostream& FunctionalPolynomial::toLatex(ostream& os, bool paren) const
{
  os << toXML();
  return os;
}

XMLObject FunctionalPolynomial::toXML() const
{
  XMLObject rtn("Polynomial");
  for (Map<MultiSet<OrderedPair<int, MultiIndex> >, 
         RefCountPtr<ScalarExpr> >::const_iterator 
         i = coeffs_.begin(); i != coeffs_.end(); i++)
    {
      const MultiSet<OrderedPair<int, MultiIndex> >& key = i->first;
      Expr e = Expr::handle(i->second);
      XMLObject term("Term");
      XMLObject coeff("Coeff");
      coeff.addChild(e.toXML());
      term.addChild(coeff);
      for (MultiSet<OrderedPair<int, MultiIndex> >::const_iterator
             j = key.begin(); j != key.end(); j++)
        {
          int funcID = j->first();
          const MultiIndex& mi = j->second();
          XMLObject f("Function");
          f.addInt("funcID", funcID);
          f.addAttribute("mi", mi.toString());
          term.addChild(f);
        }
      rtn.addChild(term);
    }
  return rtn;
}

