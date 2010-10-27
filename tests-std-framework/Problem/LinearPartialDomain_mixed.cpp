/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"
#include "SundanceCFMeshPair.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceInhomogeneousNodalDOFMap.hpp"


/** 
 * Tests logical operations on cell filters
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

CELL_PREDICATE(ATest, {return x[0] <= 0.401;})

CELL_PREDICATE(BTest, {return x[0] >= 0.399 && x[0] <= 0.601;})

CELL_PREDICATE(CTest, {return x[0] >= 0.599;})
  
//REFINE_MESH_ESTIMATE(MeshRefEst , { return 0; } , {return 1;} )

REFINE_MESH_ESTIMATE(MeshRefEst , { \
    // this is the refinement function
        if ((cellPos[0] < 0.33) && (cellPos[1] < 0.7) && (cellLevel < 1)) \
                                      return 1;\
                                 else return 0; } ,
  // this is the coarse load estimator
  { return 1; } )

MESH_DOMAIN( MeshDomain , {return true;})
  
int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      int nx = 30;
      int ny = 6;

      //MeshType meshType = new BasicSimplicialMeshType();
      //MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np, 0.0, 1.0, ny, 1, meshType);
      RefinementClass refCl = new MeshRefEst();
      MeshDomainDef meshDom = new MeshDomain();
      MeshType meshType = new HNMeshType2D();
      MeshSource mesher = new HNMesher2D(0.0, 0.0, 1.0 , 1.0 , nx , ny , meshType , refCl , meshDom );
      
      Mesh mesh = mesher.getMesh();

      Expr x = new CoordExpr(0);
      Expr dx = new Derivative(0);

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter bdry = new BoundaryCellFilter();
      CellFilter left = bdry.subset(new LeftPointTest());
      CellFilter right = bdry.subset(new RightPointTest());
      CellFilter A = interior.subset(new ATest());
      CellFilter B = interior.subset(new BTest());
      CellFilter C = interior.subset(new CTest());

    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 0);
    watchMe.setParam("discrete function evaluation", 6);
    watchMe.setParam("integration setup", 6);
    watchMe.setParam("integral transformation", 6);
    watchMe.setParam("integration", 6);
    watchMe.setParam("assembler setup", 6);
    watchMe.setParam("assembly loop", 6); 
    watchMe.setParam("matrix config",6);
    watchMe.setParam("fill", 6);
    watchMe.setParam("evaluation", 6);
    watchMe.setParam("dof map setup",6);
    watchMe.setParam("dof map access", 6); 
    //Evaluator::classVerbosity() = 6;
      
      Expr u1 = new UnknownFunction(new Lagrange(2));
      Expr u2 = new UnknownFunction(new Lagrange(1));

      Expr v1 = new TestFunction(new Lagrange(2));
      Expr v2 = new TestFunction(new Lagrange(1));

      QuadratureFamily quad = new GaussianQuadrature(2);
      Expr eqn = Integral(interior, (dx*u1)*(dx*v1), quad ) 
        + Integral(A, v2*(u2 - u1), quad ) 
        + Integral(B, v2*(u2 - 0.4), quad );
      Expr bc = EssentialBC(left, v1*u1, quad ) +
         EssentialBC(right, v1*(u1-1.0), quad );

      LinearProblem prob(mesh, eqn, bc, List(v1, v2), List(u1, u2), vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec.xml"));
#else
      ParameterXMLFileReader reader("aztec.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cout << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver); 

      // Write the field in VTK format 
      FieldWriter w = new VTKWriter("Partial_lin_Domain2d_test");
      w.addMesh(mesh);
      w.addField("u1", new ExprFieldWrapper(soln[0]));
      w.addField("u2", new ExprFieldWrapper(soln[1]));
      w.write();

      Vector<double> vec = DiscreteFunction::discFunc(soln)->getVector();

      Expr err = Integral(interior, (soln[0] - x)*(soln[0] - x) , quad )
        + Integral(A,(soln[1] - soln[0])*(soln[1] - soln[0]), quad )
        + Integral(B, (soln[1] - 0.4)*(soln[1] - 0.4), quad );

      FunctionalEvaluator errInt(mesh, err);
      double errorSq = errInt.evaluate();
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      Sundance::passFailTest(sqrt(errorSq), 1.0e-8);
      
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
