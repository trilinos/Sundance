/* @HEADER@ */
/* @HEADER@ */

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
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceStdMesh;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;



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

      int np = MPIComm::world().getNProc();

      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 2, np,
                                                         0.0, 1.0, 2, 1,
                                                         meshType);

      //      mesher.ptr()->verbosity() = VerbExtreme;

      mesher.serializeLocal();

      Mesh mesh = mesher.getMesh();

      Array<int> funcs = tuple(0);


      CellFilter interior = new MaximalCellFilter();
      CellSet interiorCells = interior.getCells(mesh);

      BasisFamily basis = new Lagrange(2);
      
      //  DOFMapBase::classVerbosity() = VerbExtreme;

      RefCountPtr<DOFMapBase> map = rcp(new HomogeneousDOFMap(mesh,
                                                              basis,
                                                              1));

      map->print(cerr);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
