#include "SundanceMeshSourceBase.hpp"

using namespace Sundance;
using namespace Sundance::StdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;


MeshSourceBase::MeshSourceBase(const MeshType& meshType,
                               const MPIComm& comm)
  : cachedMesh_(),
    hasCachedMesh_(),
    meshType_(meshType),
    comm_(comm),
    nodeAttributes_(),
    elemAttributes_()
{
}

Mesh MeshSourceBase::getMesh() const
{
  /* if we don't have a cached mesh, build one */
  if (!hasCachedMesh_)
    {
      return fillMesh();
    }
  return cachedMesh_;
}

void MeshSourceBase
::getAttributes(RefCountPtr<Array<Array<double> > >& nodeAttributes,
                RefCountPtr<Array<Array<double> > >& elemAttributes) const
{
  nodeAttributes = nodeAttributes;
  elemAttributes = elemAttributes;
}

Mesh MeshSourceBase::createMesh(int dim) const
{
  cachedMesh_ = meshType_.createEmptyMesh(dim, comm_);
  cachedMesh_.ptr()->verbosity() = verbosity();
  hasCachedMesh_ = true;
  
  return cachedMesh_;
}
