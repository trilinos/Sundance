/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBasisFamily.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

BasisFamily::BasisFamily()
  : ptr_()
{;}

BasisFamily::BasisFamily(BasisFamilyBase* ptr)
  : ptr_(ptr->getRcp())
{;}

XMLObject BasisFamily::toXML() const 
{
  return ptr_->toXML();
}

int BasisFamily::order() const 
{
  return ptr_->order();
}
