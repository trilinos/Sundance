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
{
  /* I have nonzero functional deriv of zero order, and of first order wrt 
   * myself */
  int fid = funcID();
  MultiSet<int> derivWrtMe;
  derivWrtMe.put(fid);
  addFuncIDCombo(derivWrtMe);
  addFuncIDCombo(MultiSet<int>());
}



void SymbolicFuncElement::substituteZero() const 
{
  evalPt_ = rcp(new ZeroExpr());
}

void SymbolicFuncElement
::substituteFunction(const RefCountPtr<DiscreteFuncElement>& u0) const
{
  evalPt_ = u0;
}

void SymbolicFuncElement
::findNonzeros(const EvalContext& context,
               const Set<MultiIndex>& multiIndices,
               const Set<MultiSet<int> >& activeFuncIDs,
               bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for symbolic func " 
                       << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString()
                       << " and active funcs " << activeFuncIDs);

  
  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }


  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs);

  bool isTest = (0 != dynamic_cast<const TestFuncElement*>(this));
  bool evalPtIsZero = (0 != dynamic_cast<const ZeroExpr*>(evalPt()));
  
  if (evalPtIsZero)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "eval point is a zero expr");
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tabs << "eval point is a nonzero expr");
    }

  /* Evaluate the function itself, i.e., the zeroth deriv of the function.
   * If this is a test function, or if we are doing a linear problem,
   * then we skip this step. */
  if (!regardFuncsAsConstant && !isTest && !evalPtIsZero)
    {
      if (activeFuncIDs.contains(MultiSet<int>()))
        {
          subset->addDeriv(MultipleDeriv(), VectorDeriv);
        }
      else
        {
          SUNDANCE_VERB_MEDIUM(tabs << "value of " << toString() << " not required");
        }
    }
  
  /* If this function is one of the active variables, then
   * add the deriv wrt this func to the sparsity pattern */
  MultiSet<int> myFuncID;
  myFuncID.put(funcID());
  if (activeFuncIDs.contains(myFuncID))
    {
      subset->addDeriv(new FunctionalDeriv(this, MultiIndex()),
                       ConstantDeriv);
    }
  else
    {
      SUNDANCE_VERB_MEDIUM(tabs << "deriv wrt to " << toString() << " not required");
    }
  
  const DiscreteFuncElement* df 
    = dynamic_cast<const DiscreteFuncElement*>(evalPt());
  if (df != 0)
    {
      df->findNonzeros(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant);
    }
  

  SUNDANCE_VERB_HIGH(tabs << "symbolic func " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "symbolic func " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}






