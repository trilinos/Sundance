#include "Sundance.hpp"

#include "SundanceZeroExpr.hpp"
CELL_PREDICATE(TopLeftTest, {return fabs(x[0]) < 1.0e-10 && x[1]>=0.5;}) 
CELL_PREDICATE(BottomLeftTest, {return fabs(x[0]) < 1.0e-10 && x[1]<=0.5;}) 
CELL_PREDICATE(TopRightTest, {return fabs(x[0]-1.0) < 1.0e-10 && x[1]>=0.5;}) 
CELL_PREDICATE(BottomRightTest, {return fabs(x[0]-1.0) < 1.0e-10 && x[1]<=0.5;}) 

CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;}) 

CELL_PREDICATE(InterfaceTest, {return fabs(x[1]-0.5) < 1.0e-10;}) 

CELL_PREDICATE(BottomZoneTest, {return x[1] <= 0.5;}) 
CELL_PREDICATE(TopZoneTest, {return x[1] >= 0.5;}) 


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEST_FOR_EXCEPT(npx < 1);
      TEST_FOR_EXCEPT(npy < 1);
      TEST_FOR_EXCEPT(npx * npy != np);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 32;
      int ny = 32;
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
        0.0,  1.0, ny, npy, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter(); 
      CellFilter faces = new DimensionalCellFilter(1);
      
      CellFilter topZone = interior.subset(new TopZoneTest());
      CellFilter bottomZone = interior.subset(new BottomZoneTest());
      CellFilter interface = faces.subset(new InterfaceTest());
      CellFilter topLeft = faces.subset(new TopLeftTest());
      CellFilter bottomLeft = faces.subset(new BottomLeftTest());
      CellFilter topRight = faces.subset(new TopRightTest());
      CellFilter bottomRight = faces.subset(new BottomRightTest());

      CellFilter top = faces.subset(new TopPointTest());
      CellFilter bottom = faces.subset(new BottomPointTest());

      BasisFamily basis = new Lagrange(1);
      Array<CellFilter> domains(2);
      domains[0] = top;
      domains[1] = bottom;
      DiscreteSpace discSpace(mesh, List(basis, basis), domains, vecType);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      L2Projector proj(discSpace, List(x, y));
      Expr u = proj.project();

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("twoZoneDiscFunc");
      w.addMesh(mesh);
      w.addField("x", new ExprFieldWrapper(u[0]));
      w.addField("y", new ExprFieldWrapper(u[1]));
      w.write();

      double errorSq = 0.0;
      double tol = 1.0e-6;
      Sundance::passFailTest(::sqrt(errorSq), tol);

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  Sundance::finalize();
  return Sundance::testStatus(); 
}
