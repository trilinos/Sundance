/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDimensionalCellFilter.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DimensionalCellFilter::DimensionalCellFilter(int dim)
  : CellFilterBase(), dim_(dim)
{;}

XMLObject DimensionalCellFilter::toXML() const 
{
  XMLObject rtn("DimensionalCellFilter");
  rtn.addAttribute("dim", Teuchos::toString(dim_));
  return rtn;
}

bool DimensionalCellFilter::lessThan(const CellFilterStub* other) const
{
  TEST_FOR_EXCEPTION(dynamic_cast<const DimensionalCellFilter*>(other) == 0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to DimensionalCellFilter::lessThan() should be "
                     "a DimensionalCellFilter pointer.");

  return dim_ < dynamic_cast<const DimensionalCellFilter*>(other)->dim_;
}

CellSet DimensionalCellFilter::internalGetCells(const Mesh& mesh) const
{
  return new ImplicitCellSet(mesh, dim_, 
                             mesh.cellType(dim_));
}
