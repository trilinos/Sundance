/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExplicitCellSet.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

ExplicitCellSet::ExplicitCellSet(const Mesh& mesh, int cellDim,
                                 const CellType& cellType)
  : CellSetBase(mesh, cellDim, cellType),
    cells_()
{;}

CellIterator ExplicitCellSet::begin() const
{
  return CellIterator(&cells_, CellIterator::Begin);
}

CellIterator ExplicitCellSet::end() const
{
  return CellIterator(&cells_, CellIterator::End);
}

void ExplicitCellSet::print(ostream& os) const 
{
  os << "ExplicitCellSet[cells=" << cells_ << "]";
}
