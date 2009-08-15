#include "SundanceMeshTransformationBase.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;

using namespace Teuchos;
using namespace SundanceUtils;

Mesh MeshTransformationBase::createMesh(int dim, const MPIComm& comm) const 
{
  return meshType_.createEmptyMesh(dim, comm);
}


