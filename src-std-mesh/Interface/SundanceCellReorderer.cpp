/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellReorderer.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;

RefCountPtr<CellReordererImplemBase>
CellReorderer::createInstance(const MeshBase* mesh) const
{
  return ptr()->createInstance(mesh);
}
