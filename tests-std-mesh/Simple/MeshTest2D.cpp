#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceVTKWriter.hpp"


using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;


static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}



int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 2, 2,
                                                         0.0, 1.0, 2, 1,
                                                         meshType);

      mesher.ptr()->verbosity() = VerbExtreme;
      mesher.serializeLocal();

      Mesh mesh = mesher.getMesh();

      FieldWriter w = new VerboseFieldWriter();
      FieldWriter w1 = new VerboseFieldWriter("test2d");
      FieldWriter w2 = new TriangleWriter("test2d");
      FieldWriter w3 = new VTKWriter("test2d");

      mesh.verbosity() = VerbExtreme;

      w.addMesh(mesh);
      w1.addMesh(mesh);
      w2.addMesh(mesh);
      w3.addMesh(mesh);

      w.write();
      w1.write();
      w2.write();
      w3.write();
      
    }
	catch(exception& e)
		{
      cerr << "Detected exception: " << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}