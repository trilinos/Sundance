/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMaximalCellFilter.hpp"
#include "SundanceImplicitCellSet.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

MaximalCellFilter::MaximalCellFilter()
  : CellFilterBase()
{;}

CellSet MaximalCellFilter::internalGetCells(const Mesh& mesh) const
{
  return new ImplicitCellSet(mesh, mesh.spatialDim(), 
                             mesh.cellType(mesh.spatialDim()));
}

int MaximalCellFilter::dimension(const Mesh& mesh) const 
{
  return mesh.spatialDim();
}


bool MaximalCellFilter::lessThan(const CellFilterStub* /* other */) const
{
  /* Maximal cell sets always win a comparison */
  return false;
}



