/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellPredicate.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

bool CellPredicate::operator<(const CellPredicate& other) const
{
  cerr << "comparing cell predicates \nme=" << toXML()
       << endl << "you=" << other.toXML() << endl;
  if (ptr()->typeName() < other.ptr()->typeName()) return true;
  if (ptr()->typeName() > other.ptr()->typeName()) return false;

  return ptr()->lessThan(other.ptr().get());
}


