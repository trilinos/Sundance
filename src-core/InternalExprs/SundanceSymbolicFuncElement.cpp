/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


SymbolicFuncElement::SymbolicFuncElement(const string& name,
                                         const string& suffix,
                                         int myIndex)
	: EvaluatableExpr(), FuncElementBase(name, suffix),
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

void SymbolicFuncElement::findNonzeros(const EvalContext& context,
                                       const Set<MultiIndex>& multiIndices,
                                       const Set<MultiSet<int> >& activeFuncIDs,
                                       const Set<int>& allFuncIDs,
                                       bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for symbolic func " 
                       << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString()
                       << " and active funcs " << activeFuncIDs);

  
  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }


  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);

  bool isTest = (0 != dynamic_cast<const TestFuncElement*>(this));
  bool evalPtIsZero = (0 != dynamic_cast<const ZeroExpr*>(evalPt()));
  
  SUNDANCE_VERB_MEDIUM(tabs << "eval point is zero = " << evalPtIsZero);

  /* Evaluate the function itself, i.e., the zeroth deriv of the function.
   * If this is a test function, or if we are doing a linear problem,
   * then we skip this step. */
  if (!regardFuncsAsConstant && !isTest && !evalPtIsZero)
    {
      subset->addDeriv(MultipleDeriv(), VectorDeriv);
    }
  
  /* If this function is one of the active variables, then
   * add the deriv wrt this func to the sparsity pattern */
  if (isInActiveSet(activeFuncIDs))
    {
      subset->addDeriv(new FunctionalDeriv(this, MultiIndex()),
                       ConstantDeriv);
    }

  const DiscreteFuncElement* df 
    = dynamic_cast<const DiscreteFuncElement*>(evalPt());
  if (df != 0)
    {
      df->findNonzeros(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant);
    }
  
  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  allFuncIDs, regardFuncsAsConstant);
}






