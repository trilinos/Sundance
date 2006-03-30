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

#include "SundanceDiffOp.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

DiffOp::DiffOp(const MultiIndex& op, const RefCountPtr<ScalarExpr>& arg)
  : UnaryExpr(arg), mi_(op), myCoordDeriv_(), requiredFunctions_(),
    ignoreFuncTerms_(false)
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "forming DiffOp " << toString());
  typedef Set<int>::const_iterator setIter;
  myCoordDeriv_ = new CoordDeriv(mi_.firstOrderDirection());
  
  if (isEvaluatable(arg.get()))
    {
      for (int d=0; d<MultiIndex::maxDim(); d++) 
        {
          if (op[d] == 0) continue;
          /* if the arg is constant wrt the diff op, 
           * it remains constant (zero) */ 
          if (evaluatableArg()->orderOfSpatialDependency(d) == 0 ) continue;
          /* if the arg is nonpolynomial wrt the diff op's direction, 
           * it remains
           * nonpolynomial */
          if (evaluatableArg()->orderOfSpatialDependency(d) < 0 ) 
            {
              setOrderOfDependency(d, -1);
            }
          else
            {
              /* otherwise adjust the order for differentiation */
              setOrderOfDependency(d, max(evaluatableArg()->orderOfSpatialDependency(d) - op[d], 0));
            }
        }

      setFuncIDSet(evaluatableArg()->funcIDSet());
    }
}

ostream& DiffOp::toText(ostream& os, bool /* paren */) const 
{
  string miStr = CoordExpr::coordName(mi_.firstOrderDirection(), "");
	os << "D[" << arg().toString() << ", " << miStr << "]";
	return os;
}

ostream& DiffOp::toLatex(ostream& os, bool /* paren */) const 
{
	os << "D^{" << mi_.toString() << "}" << arg().toLatex();
	return os;
}

XMLObject DiffOp::toXML() const 
{
	XMLObject rtn("DiffOp");
	rtn.addAttribute("m", mi_.toString());
  rtn.addChild(arg().toXML());

	return rtn;
}


Evaluator* DiffOp::createEvaluator(const EvaluatableExpr* expr,
                                   const EvalContext& context) const
{
  return new DiffOpEvaluator(dynamic_cast<const DiffOp*>(expr), context);
}

Set<MultiIndex> DiffOp
::argMultiIndices(const Set<MultiIndex>& multiIndices) const
{
  Set<MultiIndex> rtn;
  for (Set<MultiIndex>::const_iterator iter=multiIndices.begin();
       iter != multiIndices.end(); iter++)
    {
      rtn.put((*iter) + mi_);
    }

  return rtn;
}

Set<MultiSet<int> > DiffOp
::argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs,
                 int maxOrder) const
{
  Tabs tabs;
  Set<MultiSet<int> > rtn;// = activeFuncIDs;
  //  if (activeFuncIDs.contains(MultiSet<int>())) rtn.put(MultiSet<int>());


  SUNDANCE_VERB_MEDIUM(tabs << "DiffOp arg dependencies are " 
                       << evaluatableArg()->funcDependencies().toString());

  Set<int> allActiveFuncs;
  for (Set<MultiSet<int> >::const_iterator 
         i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& d = *i;
      for (MultiSet<int>::const_iterator j=d.begin(); j != d.end(); j++)
        {
          allActiveFuncs.put(*j);
        }
    }

  for (Set<MultiSet<int> >::const_iterator 
         i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& d = *i;
      SUNDANCE_VERB_MEDIUM(tabs << "DiffOp deriv from outside is " << d.toString());
      if (d.size() >= maxFuncDiffOrder()) continue;
      for (Set<int>::const_iterator 
             j=evaluatableArg()->funcDependencies().begin(); 
           j != evaluatableArg()->funcDependencies().end(); j++)
        {
          Tabs tab1;
          MultiSet<int> newDeriv = d;
          if (d.size() > 0 && !allActiveFuncs.contains(*j)) continue;
          SUNDANCE_VERB_MEDIUM(tab1 << "DiffOp internal dependency is " << *j);
          newDeriv.put(*j);
          SUNDANCE_VERB_MEDIUM(tab1 << "DiffOp created new arg deriv " << newDeriv.toString());
          rtn.put(newDeriv);
        }
    }
  SUNDANCE_VERB_MEDIUM(tabs << "DiffOp found arg active funcs " << rtn.toString());
  
  return evaluatableArg()->filterActiveFuncs(rtn);
}


void DiffOp::findNonzeros(const EvalContext& context,
                          const Set<MultiIndex>& multiIndices,
                          const Set<MultiSet<int> >& inputActiveFuncIDs,
                          bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for diff op " << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString());

  Set<MultiSet<int> > activeFuncIDs = filterActiveFuncs(inputActiveFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  addActiveFuncs(context, activeFuncIDs);
  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs, false);

  ignoreFuncTerms_ = regardFuncsAsConstant;
  if (ignoreFuncTerms_)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "evaluating symbolic functions at zero");
    }



  /* Figure out the sparsity pattern for the argument.  */
  Set<MultiIndex> argMI = argMultiIndices(multiIndices);
  
  SUNDANCE_VERB_MEDIUM(tabs << "DiffOp arg multi index set is " << endl << argMI);

  

  Set<MultiSet<int> > argFuncs = argActiveFuncs(activeFuncIDs, -1);

  SUNDANCE_VERB_MEDIUM(tabs << "DiffOp arg active func set is " << endl << argFuncs);

  evaluatableArg()->findNonzeros(context, argMI,
                                 argFuncs,
                                 regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> argSparsity
    = evaluatableArg()->sparsitySubset(context, argMI, argFuncs, true);

  SUNDANCE_VERB_MEDIUM(tabs << "DiffOp arg sparsity subset is " 
                       << endl << *argSparsity);


  for (int i=0; i<argSparsity->numDerivs(); i++)
    {
      Tabs tab1;

      const MultipleDeriv& md = argSparsity->deriv(i);

      if (md.order()==0) continue;
      
      SUNDANCE_VERB_MEDIUM(tab1 << "finding the effect of the argument's "
                           "nonzero derivative " << md);
      // cerr << "inverting " << md << endl;
      


      Map<MultipleDeriv, DerivState> isolatedTerms;
      Map<MultipleDeriv, Deriv> funcTerms;
      getResultDerivs(argSparsity->deriv(i), 
                      argSparsity->state(i),
                      isolatedTerms,
                      funcTerms);

      SUNDANCE_VERB_MEDIUM(tab1 << "monomials = " 
                           << isolatedTerms);


      for (Map<MultipleDeriv, DerivState>::const_iterator 
             iter=isolatedTerms.begin(); iter != isolatedTerms.end(); iter++)
        {
          subset->addDeriv(iter->first, iter->second);
        }

      if (!ignoreFuncTerms())
        {
          SUNDANCE_VERB_MEDIUM(tab1 << "func terms = " << funcTerms);
          for (Map<MultipleDeriv, Deriv>::const_iterator 
                 iter=funcTerms.begin(); iter != funcTerms.end(); iter++)
            {
              subset->addDeriv(iter->first, VectorDeriv);
            }
        }
    }


  SUNDANCE_VERB_HIGH(tabs << "diff op " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "diff op " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
  //  cerr << "my sparsity: " << *subset << endl;
}

void DiffOp::getResultDerivs(const MultipleDeriv& argDeriv,
                             const DerivState& argDerivState,
                             Map<MultipleDeriv, DerivState>& isolatedTerms,
                             Map<MultipleDeriv, Deriv>& funcTerms) const
{

  /* List the functional (not coord) single derivs in the 
   * specified multiple derivative of the argument */
  Array<Deriv> argSingleDerivs;
  for (MultipleDeriv::const_iterator j=argDeriv.begin(); 
       j != argDeriv.end(); j++)
    {
      const Deriv& d = *j;
      if (d.isCoordDeriv()) continue;
      argSingleDerivs.append(d);
    }

  TEST_FOR_EXCEPTION(argDeriv.order() > 3, RuntimeError,
                     "DiffOp::getResultDerivs() cannot handle derivatives "
                     "of order > 3");

  
  

  if (argDeriv.order()==1)
    {
      Tabs tab1;
      /* Map first-order derivatives of the argument to derivatives
       * of the diff op */

      if (argSingleDerivs.size()==0) 
        {
          /* If the argument derivative is spatial, it maps to the application
           * of the diff op to the argument, i.e., the 
           * zeroth functional derivative of the diff op expression. */
          isolatedTerms.put(MultipleDeriv(), argDerivState);
        }
      else
        {
          /* The chain rule for spatial derivatives produces a term
           * involving a first-order functional deriv times a higher
           * spatial derivative of the argument of that functional deriv.
           * Thus a first-order functional deriv of the argument
           * maps back to a zeroth order functional deriv of the diff op,
           * with a coefficient function.  */
          const FuncElementBase* f = argSingleDerivs[0].funcDeriv()->func();
          const SymbolicFuncElement* s = dynamic_cast<const SymbolicFuncElement*>(f);
          //          if (t == 0 && !ignoreFuncTerms())
          if (!s->evalPtIsZero())
            {
              Deriv d = argSingleDerivs[0].funcDeriv()->derivWrtMultiIndex(mi_);
              funcTerms.put(MultipleDeriv(), d);
            }
          
          /* The only other way a first functional deriv will appear is 
           * through first functional differentiation of the spatial
           * chain rule. We need to back out the variable of differentiation
           * that produced this term. It is possible that no such variable
           * exists */
          
          Deriv dBack;
          if (canBackOutDeriv(argSingleDerivs[0], mi_, dBack))
            {
              MultipleDeriv mdBack;
              mdBack.put(dBack);
              isolatedTerms.put(mdBack, argDerivState);
            }
        }
    }
  else if (argDeriv.order()==2)
    {
      /* Map second-order derivatives of the argument to derivatives
       * of the diff op */
      
      if (argSingleDerivs.size()==1)
        {
          /* A first-order functional derivative here means we have a mixed
           * spatial-functional second derivative. That corresponds to 
           * the "commuting" term in a first functional deriv of the diff op.
           * The variable of differentiation here is the single deriv. */
          MultipleDeriv md1;
          md1.put(argSingleDerivs[0]);
          isolatedTerms.put(md1, argDerivState);
        }
      else 
        {
          /* Doing functional differentiation on a chain-rule expansion
           * of a spatial derivative produces terms
           * involving second-order functional derivs times a higher
           * spatial derivative of the function that is the argument of
           * functional differentiation. */
          MultipleDeriv md1;
          md1.put(argSingleDerivs[0]);
          MultipleDeriv md2;
          md2.put(argSingleDerivs[1]);

          const FuncElementBase* f0 = argSingleDerivs[0].funcDeriv()->func();
          const SymbolicFuncElement* s0 = dynamic_cast<const SymbolicFuncElement*>(f0);

          const FuncElementBase* f1 = argSingleDerivs[1].funcDeriv()->func();
          const SymbolicFuncElement* s1 = dynamic_cast<const SymbolicFuncElement*>(f1);

          //          if (t0==0 && !ignoreFuncTerms())
          if (!s0->evalPtIsZero())
            {
              Deriv d1 = argSingleDerivs[0].funcDeriv()->derivWrtMultiIndex(mi_);
              funcTerms.put(md2, d1);
            }

          //         if (t1==0 && !ignoreFuncTerms())
          if (!s1->evalPtIsZero())
            {
              Deriv d2 = argSingleDerivs[1].funcDeriv()->derivWrtMultiIndex(mi_);

              funcTerms.put(md1, d2);
            }

          
          
          Deriv dBack;
          if (canBackOutDeriv(argSingleDerivs[0], mi_, dBack))
            {
              MultipleDeriv md1;
              md1.put(argSingleDerivs[1]);
              md1.put(dBack);
              isolatedTerms.put(md1, argDerivState);
            }
          if (canBackOutDeriv(argSingleDerivs[1], mi_, dBack))
            {
              MultipleDeriv md1;
              md1.put(argSingleDerivs[0]);
              md1.put(dBack);
              isolatedTerms.put(md1, argDerivState);
            }
        }
    }
  else if (argDeriv.order()==3)
    {
      if (argSingleDerivs.size()==2)
        {
          /* A mixed second-order functional, first-order spatial deriv here is
           * the commuting term in the chain rule. */
          MultipleDeriv md1;
          md1.put(argSingleDerivs[0]);
          md1.put(argSingleDerivs[1]);
          isolatedTerms.put(md1, argDerivState);
        }
      else
        {
          /* a third-order functional 
           * deriv appears in the chain rule expansion of a second deriv */
          MultipleDeriv md12;
          MultipleDeriv md13;
          MultipleDeriv md23;

          md12.put(argSingleDerivs[0]);
          md12.put(argSingleDerivs[1]);

          md13.put(argSingleDerivs[0]);
          md13.put(argSingleDerivs[2]);

          md23.put(argSingleDerivs[1]);
          md23.put(argSingleDerivs[2]);

          

          Deriv d1 = argSingleDerivs[0].funcDeriv()->derivWrtMultiIndex(mi_);
          Deriv d2 = argSingleDerivs[1].funcDeriv()->derivWrtMultiIndex(mi_);
          Deriv d3 = argSingleDerivs[2].funcDeriv()->derivWrtMultiIndex(mi_);

          funcTerms.put(md12, d3);
          funcTerms.put(md13, d2);
          funcTerms.put(md23, d1);
        }
    }
}


bool DiffOp::canBackOutDeriv(const Deriv& d, const MultiIndex& x, 
                             Deriv& rtnDeriv) const
{
  TEST_FOR_EXCEPTION(d.isCoordDeriv(), InternalError,
                     "DiffOp::canBackOutDeriv should not be called for "
                     "spatial derivative");

  MultiIndex alpha = d.funcDeriv()->multiIndex();
  MultiIndex miNew;
  for (int i=0; i<MultiIndex::maxDim(); i++)
    {
      miNew[i] = alpha[i] + x[i];
      if (miNew[i] < 0) return false;
    }
  rtnDeriv = new FunctionalDeriv(d.funcDeriv()->func(), miNew);
  return true;
}



RefCountPtr<Array<Set<MultipleDeriv> > > 
DiffOp::internalDetermineR(const EvalContext& context,
                           const Array<Set<MultipleDeriv> >& RInput) const
{
  const Deriv& x = myCoordDeriv();

  /* record R for self */
  RefCountPtr<Array<Set<MultipleDeriv> > > rtn 
    = EvaluatableExpr::internalDetermineR(context, RInput);

   
  /* work out RInput for the argument */
  Array<Set<MultipleDeriv> > RArg(RInput.size()+1);

  for (unsigned int order=0; order<RArg.size(); order++)
    {
      const Set<MultipleDeriv>& Warg = evaluatableArg()->findW(order, context);      
      const Set<MultipleDeriv>& R = RInput[order];
      Set<MultipleDeriv> RMinus;
      if (order > 0) RMinus = RInput[order-1];
      Set<MultipleDeriv>& tmp = RArg[order];
      
      Set<MultipleDeriv>::const_iterator i;
      for (i=Warg.begin(); i!=Warg.end(); i++)
        {
          if (!(i->contains(x))) continue;
          MultipleDeriv xRemoved = i->factorOutDeriv(x);
          if (RMinus.contains(xRemoved)) tmp.put(*i);
        }
  
      for (i=Warg.begin(); i!=Warg.end(); i++)
        {
          const MultipleDeriv& md = *i;
          for (MultipleDeriv::const_iterator j=md.begin(); j!=md.end(); j++)
            {
              const Deriv& d = *j;
              /* coord derivs are handled elsewhere, so skip any coord derivs */
              if (!d.isFunctionalDeriv()) continue;
              /* see if the multiple deriv can be "factored" into a single deriv
               * corresponding to one of the nonzero functions, plus a remaining 
               * multiple deriv that is a member of RMinus */
              MultipleDeriv remainder = md.factorOutDeriv(d);
              const FunctionalDeriv* fd = d.funcDeriv();
              const FuncElementBase* f = fd->func();
              const SymbolicFuncElement* s = dynamic_cast<const SymbolicFuncElement*>(f);
              TEST_FOR_EXCEPTION(s==0, RuntimeError, "non-symbolic function in functional "
                                 "derivative " << d);
          
              if (!s->evalPtIsZero() && RMinus.contains(remainder))
                {
                  tmp.put(md);
                  /* If we've added this multiple deriv to the set, there's no need to
                   * test the other permutations so we can skip the rest of the loop. */
                  break;
                }
            }
        }

      for (i=Warg.begin(); i!=Warg.end(); i++)
        {
          Tabs tab1;
          const MultipleDeriv& md = *i;

          for (MultipleDeriv::const_iterator j=md.begin(); j!=md.end(); j++)
            {
              Tabs tab2;
              const Deriv& d = *j;
              if (!d.isFunctionalDeriv()) continue;
              MultipleDeriv remainder = md.factorOutDeriv(d);
              remainder.put(d.funcDeriv()->derivWrtMultiIndex(mi()));
              if (R.contains(remainder)) tmp.put(remainder);
            }
        }
    }

  /* Propagate the computation of R into the argument */
  evaluatableArg()->determineR(context, RArg);

  
  /* all done */  
  return rtn;
}


Set<MultipleDeriv> DiffOp::internalFindW(int order, const EvalContext& context) const
{
  Tabs tabs;
  const Deriv& x = myCoordDeriv();
  const Set<MultipleDeriv>& WArg = evaluatableArg()->findW(order, context);
  const Set<MultipleDeriv>& WArgPlus = evaluatableArg()->findW(order+1, context);
  Set<MultipleDeriv> rtn;

  Set<MultipleDeriv>::const_iterator i;
  
  /* */
  cout << tabs << "W(order) for arg = " << WArg << endl;
  cout << tabs << "W(order+1) for arg = " << WArgPlus << endl;
  for (i=WArgPlus.begin(); i!=WArgPlus.end(); i++)
    {
      Tabs tab1;
      /* look for a nonzero mixed spatial/functional deriv */
      const MultipleDeriv& md = *i;
      cout << tab1 << "checking arg deriv = " << md << endl;
      if (md.contains(x)) 
        {
          Tabs tab2;
          MultipleDeriv remainder = md.factorOutDeriv(x);
          cout << tab2 << "adding = " << remainder << endl;
          rtn.put(md.factorOutDeriv(x));
        }

      /* look for a nonzero functional deriv times a nonzero function */ 
      for (MultipleDeriv::const_iterator j=md.begin(); j!=md.end(); j++)
        {
          Tabs tab2;
          const Deriv& d = *j;
          if (!d.isFunctionalDeriv()) continue;
          MultipleDeriv remainder = md.factorOutDeriv(d);
          const FunctionalDeriv* fd = d.funcDeriv();
          const FuncElementBase* f = fd->func();
          const SymbolicFuncElement* s = dynamic_cast<const SymbolicFuncElement*>(f);
          TEST_FOR_EXCEPTION(s==0, RuntimeError, "non-symbolic function in functional "
                             "derivative " << d);
          if (!s->evalPtIsZero())
            {
              cout << tab2 << "adding = " << remainder << endl;
              rtn.put(remainder);
            }
        }
    }

  /* */
  cout << tabs << "doing delta terms" << endl;
  for (i=WArg.begin(); i!=WArg.end(); i++)
    {
      Tabs tab1;
      const MultipleDeriv& md = *i;
      cout << tab1 << "checking arg deriv = " << md << endl;

      for (MultipleDeriv::const_iterator j=md.begin(); j!=md.end(); j++)
        {
          Tabs tab2;
          const Deriv& d = *j;
          if (!d.isFunctionalDeriv()) continue;
          MultipleDeriv remainder = md.factorOutDeriv(d);
          remainder.put(d.funcDeriv()->derivWrtMultiIndex(mi()));
          cout << tab2 << "adding = " << remainder << endl;
          rtn.put(remainder);
        }
    }

  return rtn;
}


