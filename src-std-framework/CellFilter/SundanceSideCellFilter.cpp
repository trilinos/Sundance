/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSideCellFilter.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

SideCellFilter::SideCellFilter()
  : CellFilterBase()
{;}

CellSet SideCellFilter::internalGetCells(const Mesh& mesh) const
{
  return new ImplicitCellSet(mesh, mesh.spatialDim()-1, 
                             mesh.cellType(mesh.spatialDim()-1));
}

int SideCellFilter::dimension(const Mesh& mesh) const 
{
  return mesh.spatialDim()-1;
}


bool SideCellFilter::lessThan(const CellFilterStub* other) const
{
  const SideCellFilter* S 
    = dynamic_cast<const SideCellFilter*>(other);

  TEST_FOR_EXCEPTION(S==0,
                     InternalError,
                     "argument " << other->toXML() 
                     << " to SideCellFilter::lessThan() should be "
                     "a SideCellFilter pointer.");
}



