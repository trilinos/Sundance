#include "SundanceMeshSourceBase.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;


MeshSourceBase::MeshSourceBase(const MeshType& meshType,
                               const MPIComm& comm)
  : cachedMesh_(),
    hasCachedMesh_(),
    meshType_(meshType),
    comm_(comm),
    nodeAttributes_(rcp(new Array<Array<double> >())),
    elemAttributes_(rcp(new Array<Array<double> >()))
{
}

Mesh MeshSourceBase::getMesh() const
{
  Tabs tabs;
  
  /* if we don't have a cached mesh, build one */
  if (!hasCachedMesh_)
    {
      Mesh rtn =  fillMesh();
      cerr << tabs << "got a mesh with " << rtn.numCells(0)
                        << " nodes and " << rtn.numCells(rtn.spatialDim())
           << " maximal cells" << endl;
      return rtn;
    }
  return cachedMesh_;
}

void MeshSourceBase
::getAttributes(RefCountPtr<Array<Array<double> > >& nodeAttributes,
                RefCountPtr<Array<Array<double> > >& elemAttributes) const
{
  nodeAttributes = nodeAttributes_;
  elemAttributes = elemAttributes_;
}

Mesh MeshSourceBase::createMesh(int dim) const
{
  cachedMesh_ = meshType_.createEmptyMesh(dim, comm_);
  cachedMesh_.ptr()->verbosity() = verbosity();
  hasCachedMesh_ = true;
  
  return cachedMesh_;
}
