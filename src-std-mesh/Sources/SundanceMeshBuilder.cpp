#include "SundanceMeshBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceTriangleMeshReader.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"



using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;



Mesh MeshBuilder::createMesh(const ParameterList& params)
{
  TEST_FOR_EXCEPTION(!params.isParameter("type"), RuntimeError,
                     "field name 'type' expected but not found in MeshBuilder "
                     "input parameter list: " << params);

  string type = params.get<string>("type");

  MeshSource mesher;

  if (type=="Rectangle")
    {
      mesher = new PartitionedRectangleMesher(params);
    }
  else if (type=="Line")
    {
      mesher = new PartitionedLineMesher(params);
    }
  else if (type=="Exodus")
    {
      mesher = new ExodusNetCDFMeshReader(params);
    }
  else if (type=="Triangle")
    {
      mesher = new ExodusNetCDFMeshReader(params);
    }

  TEST_FOR_EXCEPTION(mesher.ptr().get()==0, RuntimeError,
                     "null mesh source in MeshBuilder::createMesh()");

  return mesher.getMesh();
}