/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCoordExpr.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CoordExpr::CoordExpr(int dir, const string& name)
  : LeafExpr(), FuncElementBase(coordName(dir, name), ""), 
    dir_(dir)
{
  setOrderOfDependency(dir, 1);
}

XMLObject CoordExpr::toXML() const 
{
  XMLObject rtn("CoordExpr");
  rtn.addAttribute("dir", Teuchos::toString(dir_));
  rtn.addAttribute("name", name());
  return rtn;
}

string CoordExpr::coordName(int dir, const string& name)
{
  if (name.length() > 0) return name;
  switch(dir)
    {
    case 0:
      return "x";
    case 1:
      return "y";
    case 2:
      return "z";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "CoordExpr::coordName direction out of range [0,2]");
      return "error";
    }
}



void CoordExpr::findNonzeros(const EvalContext& context,
                             const Set<MultiIndex>& multiIndices,
                                const Set<MultiSet<int> >& activeFuncIDs,
                                const Set<int>& allFuncIDs,
                             bool regardFuncsAsConstant) const
{
  Tabs tabs;
  SUNDANCE_VERB_MEDIUM(tabs << "finding nonzeros for coord func " << toString()
                       << " subject to multi index set " 
                       << multiIndices.toString());

  if (nonzerosAreKnown(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant))
    {
      SUNDANCE_VERB_MEDIUM(tabs << "...reusing previously computed data");
      return;
    }

  RefCountPtr<SparsitySubset> subset = sparsitySubset(context, multiIndices);

  MultiIndex myDirection;
  myDirection[dir_] = 1;
  if (multiIndices.contains(myDirection)) 
    {
      subset->addDeriv(new CoordDeriv(dir_), ConstantDeriv);
    }
  MultiIndex empty;
  if (multiIndices.contains(empty))
    {
      subset->addDeriv(MultipleDeriv(), VectorDeriv);
    }

  addKnownNonzero(context, multiIndices, activeFuncIDs,
                       allFuncIDs, regardFuncsAsConstant);
}




