#include "SundanceMeshTransformationBase.hpp"

using namespace Sundance;
using namespace Sundance::StdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;

Mesh MeshTransformationBase::createMesh(int dim, const MPIComm& comm) const 
{
  return meshType_.createEmptyMesh(dim, comm);
}


