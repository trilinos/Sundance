#include "SundanceMeshTransformation.hpp"

using namespace Sundance;
using namespace Sundance::StdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;



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
  if (rtn.spatialDim() > 1) rtn.assignIntermediateCellOwners(1);
  if (rtn.spatialDim() > 2) rtn.assignIntermediateCellOwners(2);
  return rtn;
}
