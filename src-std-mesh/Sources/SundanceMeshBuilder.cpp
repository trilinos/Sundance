#include "SundanceMeshBuilder.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"

#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceTriangleMeshReader.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"



using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;



Mesh MeshBuilder::createMesh(const ParameterList& params)
{
  TEST_FOR_EXCEPTION(!params.isParameter("type"), std::runtime_error,
                     "field name 'type' expected but not found in MeshBuilder "
                     "input parameter list: " << params);

  std::string type = params.get<string>("type");

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
      mesher = new TriangleMeshReader(params);
    }

  TEST_FOR_EXCEPTION(mesher.ptr().get()==0, std::runtime_error,
                     "null mesh source in MeshBuilder::createMesh()");

  return mesher.getMesh();
}
