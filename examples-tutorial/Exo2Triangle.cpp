#include "Sundance.hpp"

/** 
 * Exo2Triangle converts a NCDF mesh to a Triangle format mesh
 */


int main(int argc, void** argv)
{
  
  try
		{
      string inputFile;
      string outputFile;
      int indexOffset=0;

      Sundance::setOption("i", inputFile, "Input mesh file");
      Sundance::setOption("o", outputFile, "Output mesh file");
      Sundance::setOption("offset", indexOffset, "Index offset");

      Sundance::init(&argc, &argv);

      MeshType meshType = new BasicSimplicialMeshType();

      /* Read the input mesh */
      MeshSource mesher 
        = new ExodusNetCDFMeshReader(inputFile, meshType);
      Mesh mesh = mesher.getMesh();


      FieldWriter w = new TriangleWriter(outputFile, indexOffset);
      w.addMesh(mesh);
      w.write();


    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  TimeMonitor::summarize();
  Sundance::finalize();
}
