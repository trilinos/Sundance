/* @HEADER@ */
/* @HEADER@ */

#include "SundanceLabeledSubsetCellFilter.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

LabeledSubsetCellFilter::LabeledSubsetCellFilter(const CellFilter& superset,
                                                 const string& label)
  : SubsetCellFilter(superset, new LabelPredicate(label))
{;}


XMLObject LabeledSubsetCellFilter::toXML() const 
{
  XMLObject rtn("LabeledSubsetCellFilter");
  rtn.addAttribute("label", label_);
  return rtn;
}

bool LabeledSubsetCellFilter::lessThan(const CellFilterBase* other) const
{
  const LabeledSubsetCellFilter* L 
    = dynamic_cast<const LabeledSubsetCellFilter*>(other);

  TEST_FOR_EXCEPTION(L==0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to LabeledSubsetCellFilter::lessThan() should be "
                     "a LabeledSubsetCellFilter pointer.");

  return OrderedPair<label_, superset_> 
    < OrderedPair<L->label_, L->superset_>;
}
