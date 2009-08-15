#include "SundanceMeshTransformation.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;

using namespace Teuchos;
using namespace SundanceUtils;



MeshTransformation::MeshTransformation()
  : Handle<MeshTransformationBase>()
{}

MeshTransformation::MeshTransformation(Handleable<MeshTransformationBase>* rawPtr)
  : Handle<MeshTransformationBase>(rawPtr)
{}


MeshTransformation::MeshTransformation(const RefCountPtr<MeshTransformationBase>& smartPtr)
  : Handle<MeshTransformationBase>(smartPtr)
{}

Mesh MeshTransformation::apply(const Mesh& inputMesh) const 
{
  Mesh rtn = ptr()->apply(inputMesh);
  //if (rtn.spatialDim() > 1) rtn.assignIntermediateCellOwners(1);
  //if (rtn.spatialDim() > 2) rtn.assignIntermediateCellOwners(2);
  return rtn;
}
