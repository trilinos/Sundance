/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellDiameterExpr::CellDiameterExpr(const string& name)
  : LeafExpr(), name_(name)
{}

XMLObject CellDiameterExpr::toXML() const 
{
  XMLObject rtn("CellDiameterExpr");
  rtn.addAttribute("name", name_);
  return rtn;
}


void CellDiameterExpr::findNonzeros(const EvalContext& context,
                             const Set<MultiIndex>& multiIndices,
                             bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for cell diameter " << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString());

  if (nonzerosAreKnown(context, multiIndices, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);

  MultiIndex empty;
  if (multiIndices.contains(empty))
    {
      subset->addDeriv(MultipleDeriv(), VectorDeriv);
    }

  addKnownNonzero(context, multiIndices, regardFuncsAsConstant);
}

ostream& CellDiameterExpr::toText(ostream& os, bool paren) const
{
  os << name();
  return os;
}


ostream& CellDiameterExpr::toLatex(ostream& os, bool paren) const
{
  os << name();
  return os;
}


