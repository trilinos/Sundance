/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordDeriv.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

DiscreteFuncElement::DiscreteFuncElement(DiscreteFunctionStub* master, 
                                         const string& name,
                                         const string& suffix,
                                         int myIndex)
	: LeafExpr(), 
    FuncElementBase(name, suffix),
    master_(master),
    myIndex_(myIndex)
{}

void DiscreteFuncElement::findNonzeros(const EvalContext& context,
                                       const Set<MultiIndex>& multiIndices,
                                       const Set<MultiSet<int> >& activeFuncIDs,
                                       bool regardFuncsAsConstant) const
{

  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for discrete func " 
                       << toString() << " subject to multiindices "
                       << multiIndices);

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }


  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices, activeFuncIDs);
  
  for (Set<MultiIndex>::const_iterator 
         i=multiIndices.begin(); i != multiIndices.end(); i++)
    {
      if (i->order()==1)
        {
          subset->addDeriv(new CoordDeriv(i->firstOrderDirection()), 
                           VectorDeriv);
        }
      if (i->order()==0)
        {
          subset->addDeriv(MultipleDeriv(),
                           VectorDeriv);
        }
    }
  
  SUNDANCE_VERB_HIGH(tabs << "discrete func " + toString()
                     << ": my sparsity subset is " 
                     << endl << *subset);

  SUNDANCE_VERB_HIGH(tabs << "discrete func " + toString() 
                     << " my sparsity superset is " 
                     << endl << *sparsitySuperset(context));
  addKnownNonzero(context, multiIndices, activeFuncIDs,
                  regardFuncsAsConstant);
}

XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}

