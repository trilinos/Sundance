/* @HEADER@ */
/* @HEADER@ */

#include "SundanceImplicitCellSet.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

ImplicitCellSet::ImplicitCellSet(const Mesh& mesh, int cellDim,
                                 const CellType& cellType)
  : CellSetBase(mesh, cellDim, cellType)
{;}

CellIterator ImplicitCellSet::begin() const
{
  return CellIterator(mesh(), dimension(), CellIterator::Begin);
}

CellIterator ImplicitCellSet::end() const
{
  return CellIterator(mesh(), dimension(), CellIterator::End);
}

void ImplicitCellSet::print(ostream& os) const 
{
  os << "ImplicitCellSet[dim=" << dimension() << ", type=" << cellType() << "]";
}

