#include "SundanceMeshType.hpp"

using namespace Sundance;
using namespace Sundance::StdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;



MeshType::MeshType()
  : Handle<MeshTypeBase>()
{}

MeshType::MeshType(Handleable<MeshTypeBase>* rawPtr)
  : Handle<MeshTypeBase>(rawPtr)
{}


MeshType::MeshType(const RefCountPtr<MeshTypeBase>& smartPtr)
  : Handle<MeshTypeBase>(smartPtr)
{}

Mesh MeshType::createEmptyMesh(int dim, const MPIComm& comm) const 
{
  Mesh rtn;
  try
    {
      rtn = ptr()->createEmptyMesh(dim, comm);
    }
  catch(std::exception& e)
    {
      SUNDANCE_TRACE(e);
    }
  return rtn;
}


