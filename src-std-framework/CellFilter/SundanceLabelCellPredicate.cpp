/* @HEADER@ */
/* @HEADER@ */

#include "SundanceLabelCellPredicate.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

bool LabelCellPredicate::lessThan(const CellPredicateBase* other) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const LabelCellPredicate*>(other) == 0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to LabelCellPredicate::lessThan() should be "
                     "a LabelCellPredicate pointer.");

  return labelIndex_ < dynamic_cast<const LabelCellPredicate*>(other)->labelIndex_;
}

bool LabelCellPredicate::test(int cellLID) const 
{
  
  return mesh().label(cellDim(), cellLID) == labelIndex_;
}

XMLObject LabelCellPredicate::toXML() const 
{
  XMLObject rtn("LabelCellPredicate");
  rtn.addAttribute("label", Teuchos::toString(labelIndex_));
  return rtn;
}

