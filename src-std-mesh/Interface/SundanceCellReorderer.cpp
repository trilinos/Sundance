/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellReorderer.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance::StdMesh::Internal;
using namespace Sundance;
using namespace Teuchos;

RefCountPtr<CellReordererImplemBase>
CellReorderer::createInstance(const MeshBase* mesh) const
{
  return ptr()->createInstance(mesh);
}
