/* @HEADER@ */
/* @HEADER@ */

#include "SundanceNullCellFilter.hpp"
#include "SundanceExplicitCellSet.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

NullCellFilter::NullCellFilter()
  : CellFilterBase()
{;}

CellSet NullCellFilter::internalGetCells(const Mesh& mesh) const
{
  return new ExplicitCellSet(mesh, -1, 
                             NullCell);
}

int NullCellFilter::dimension(const Mesh& mesh) const 
{
  -1;
}

XMLObject NullCellFilter::toXML() const 
{
  XMLObject rtn("NullCellFilter");
  return rtn;
}

bool NullCellFilter::lessThan(const CellFilterBase* /* other */) const
{
  /* Maximal cell sets always lose a comparison */
  return true;
}
