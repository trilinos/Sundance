/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

SymbolicFuncElement::SymbolicFuncElement(const string& name,
                                         int myIndex)
	: FuncElementBase(name), EvaluatableExpr(), 
    evalPt_(),
    evalPtDerivSetIndices_(),
    myIndex_(myIndex)
{}



void SymbolicFuncElement::substituteZero() const 
{
  evalPt_ = rcp(new ZeroExpr());
}

void SymbolicFuncElement
::substituteFunction(const RefCountPtr<DiscreteFuncElement>& u0) const
{
  evalPt_ = u0;
}

bool SymbolicFuncElement::hasNonzeroDeriv(const MultipleDeriv& d) const
{
  TEST_FOR_EXCEPTION(evalPt_.get() == NULL, InternalError,
                     "SymbolicFuncElement::hasNonzeroDeriv() detected an unknown "
                     "function with an undefined evaluation point. "
                     "Please define an evaluation point for the function "
                     << toString());

  /* If we are evaluating a zeroth derivative, i.e., the function itself,
   * the derivative is nonzero */
  if (d.order()==0) return true;

  /* If we are evaluating a first derivative, the derivative is zero
   * unless it is a functional derivative with respect to this function */
  if (d.order()==1)
    {
      Deriv deriv = *(d.begin());
      const FunctionalDeriv* f = deriv.funcDeriv();

      return (f != 0 && funcID()==f->funcID() && 
              f->multiIndex().order()==0);
    }

  /* All higher-order functional derivatives are zero */
  return false;
}

void SymbolicFuncElement::getRoughDependencies(Set<Deriv>& funcs) const
{
  funcs.put(new FunctionalDeriv(this, MultiIndex()));
}


int SymbolicFuncElement::setupEval(const EvalContext& region,
                               const EvaluatorFactory* factory) const
{
  /* If we've been here already with this deriv set, we're done.
   * If we've been here with this deriv set, but in another region,
   * map the region to the deriv set. */
  bool derivSetIsKnown;
  if (checkForKnownRegion(region, derivSetIsKnown))
    {
      return getDerivSetIndex(region);
    }

  /* Create a new entry in our tables of deriv sets. This step creates
  * a sparsity pattern for the current deriv set. */
  int derivSetIndex = registerRegion(region, derivSetIsKnown,
                                     currentDerivSuperset(), 
                                     factory);
  
  /* set up the eval point, returning the index by which the eval point
   * refers to this deriv set.  */
  int evalPtDerivSetIndex 
    = evalPt()->setupEval(region, 
                          factory);
  
  /* store the eval point's deriv set index */
  evalPtDerivSetIndices_.append(evalPtDerivSetIndex);
  
  return derivSetIndex;
}

void SymbolicFuncElement::resetDerivSuperset() const 
{
  currentDerivSuperset() = DerivSet();
  
  evalPt()->resetDerivSuperset();
}

void SymbolicFuncElement::findDerivSuperset(const DerivSet& derivs) const 
{
  Tabs tabs;
  if (verbosity() > 1)
    {
      cerr << tabs << "finding deriv superset for symbol " << toString()
           << endl;
    }
  currentDerivSuperset().merge(derivs);

  evalPt()->findDerivSuperset(derivs);
}
