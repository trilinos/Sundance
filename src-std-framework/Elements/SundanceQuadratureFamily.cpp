/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureFamily.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

void QuadratureFamily::getPoints(const CellType& cellType, 
                                 Array<Point>& quadPoints,
                                 Array<double>& quadWeights) const 
{
  const QuadratureFamilyBase* q 
    = dynamic_cast<const QuadratureFamilyBase*>(ptr().get());
  
  TEST_FOR_EXCEPTION(q==0, InternalError, 
                     "QuadratureFamilyStub pointer" << toXML().toString() 
                     << " could not be cast to a QuadratureFamilyBase ptr");
                     
  q->getPoints(cellType, quadPoints, quadWeights);
}

XMLObject QuadratureFamily::toXML() const 
{
  return ptr()->toXML();
}

int QuadratureFamily::order() const 
{
  return ptr()->order();
}
