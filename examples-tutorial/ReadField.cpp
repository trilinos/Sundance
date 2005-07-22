#include "Sundance.hpp"

int main(int argc, void** argv)
{
  
  try
    {
      MPISession::init(&argc, &argv);
      
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create the mesh */
      int n = 4;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, n, 1,
                                                         0.0, 1.0, n, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();


      /* create an empty (zero-valued) discrete function */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");

      /* get from the discrete function a pointer to the underlying
       * vector and a pointer to the node-to-dof map. */
      Vector<double> vec = DiscreteFunction::discFunc(u0)->getVector();
      const RefCountPtr<DOFMapBase>& dofMap 
	= DiscreteFunction::discFunc(u0)->map();

      /* Read the file from whatever format. Here, we assume it's in
       * a simple text file 
       * 
       * Line 1: number of nodes
       * Lines 2 through N+1: data */
      ifstream is("../../examples-tutorial/fieldData.dat");
      int nNodes;
      is >> nNodes;

      /* check for consistency between the mesh and the data */
      TEST_FOR_EXCEPTION(mesh.numCells(0) != nNodes, RuntimeError,
			 "number of nodes in data file fieldData.dat is " 
			 << nNodes << " but number of nodes in mesh is "
			 << mesh.numCells(0));
      
      /* read the data, putting each entry into its correct place in
       * the vector as indexed by the dof number, NOT the node number.
       * Here I'm assuming one field value per node. */
      Array<int> dofs(1);
      double fVal;
      for (int i=0; i<nNodes; i++)
	{
	  /* look up the dof for the 0-th function on this node */
	  dofMap->getDOFsForCell(0, i, 0, dofs);
	  int dof = dofs[0];
	  is >> fVal;
	  vec.setElement(dof, fVal);
	}
      DiscreteFunction::discFunc(u0)->setVector(vec);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("FieldReaderTest");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(u0[0]));
      w.write();
    }
  catch(exception& e)
    {
      cerr << e.what() << endl;
    }
  Sundance::finalize();
}
