/* @HEADER@ */
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
::argActiveFuncs(const Set<MultiSet<int> >& activeFuncIDs) const
{
  Set<MultiSet<int> > rtn = activeFuncIDs;
  for (Set<MultiSet<int> >::const_iterator 
         i=activeFuncIDs.begin(); i != activeFuncIDs.end(); i++)
    {
      const MultiSet<int>& d = *i;
      if (d.size() >= maxFuncDiffOrder()) continue;
      for (Set<int>::const_iterator 
             j=funcDependencies().begin(); j != funcDependencies().end(); j++)
        {
          MultiSet<int> newDeriv = d;
          newDeriv.put(*j);
          rtn.put(newDeriv);
        }
    }
  return rtn;
}


void DiffOp::findNonzeros(const EvalContext& context,
                          const Set<MultiIndex>& multiIndices,
                          const Set<MultiSet<int> >& activeFuncIDs,
                          bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for diff op " << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString());
  SUNDANCE_VERB_MEDIUM(tabs << "active funcs are " << activeFuncIDs);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  addActiveFuncs(context, activeFuncIDs);
  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs);

  ignoreFuncTerms_ = regardFuncsAsConstant;
  if (ignoreFuncTerms_)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "evaluating symbolic functions at zero");
    }



  /* Figure out the sparsity pattern for the argument.  */
  Set<MultiIndex> argMI = argMultiIndices(multiIndices);
  
  SUNDANCE_VERB_MEDIUM(tabs << "arg multi index set is " << endl << argMI);

  

  Set<MultiSet<int> > argFuncs = argActiveFuncs(activeFuncIDs);
  evaluatableArg()->findNonzeros(context, argMI,
                                 argFuncs,
                                 regardFuncsAsConstant);

  RefCountPtr<SparsitySubset> argSparsity
    = evaluatableArg()->sparsitySubset(context, argMI, argFuncs);

  SUNDANCE_VERB_MEDIUM(tabs << "arg sparsity subset is " 
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
          const TestFuncElement* t = dynamic_cast<const TestFuncElement*>(f);
          if (t == 0 && !ignoreFuncTerms())
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
          const TestFuncElement* t0 = dynamic_cast<const TestFuncElement*>(f0);

          const FuncElementBase* f1 = argSingleDerivs[1].funcDeriv()->func();
          const TestFuncElement* t1 = dynamic_cast<const TestFuncElement*>(f1);

          if (t0==0 && !ignoreFuncTerms())
            {
              Deriv d1 = argSingleDerivs[0].funcDeriv()->derivWrtMultiIndex(mi_);
              funcTerms.put(md2, d1);
            }

          if (t1==0 && !ignoreFuncTerms())
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



