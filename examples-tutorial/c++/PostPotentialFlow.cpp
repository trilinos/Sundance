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

/** 
 * \example PostPotentialFlow.cpp
 * Solves the Laplace equation for potential flow past an elliptical 
 * post in a wind tunnel.
 */


int main(int argc, void** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("../../../examples-tutorial/meshes/post.ncdf", meshType);
      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter boundary = new BoundaryCellFilter();
      CellFilter in = boundary.labeledSubset(1);
      CellFilter out = boundary.labeledSubset(2);
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr phi = new UnknownFunction(new Lagrange(1), "u");
      Expr phiHat = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate functions */
      Expr x = new CoordExpr(0);
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);
      
      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      double L = 1.0;
      /* Define the weak form */
      Expr eqn = Integral(interior, (grad*phiHat)*(grad*phi), quad2)
        + Integral(in, phiHat*(x-phi)/L, quad2) ;


      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(out, phiHat*phi/L, quad2);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, phiHat, phi, vecType);


      /* Read the parameters for the linear solver from an XML file */
      ParameterXMLFileReader reader("../../../examples-tutorial/data/bicgstab.xml");
      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(solverParams);


      /* solve the problem */
      Expr soln = prob.solve(linSolver);


      
      // /* Project the velocity onto a discrete space so we can visualize it */
      DiscreteSpace discreteSpace(mesh, 
                                  List(new Lagrange(1), 
                                       new Lagrange(1), 
                                       new Lagrange(1)),
                                  vecType);
      L2Projector projector(discreteSpace, grad*soln);
      Expr velocity = projector.project();

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Post3d");
      w.addMesh(mesh);
      w.addField("phi", new ExprFieldWrapper(soln[0]));
      w.addField("ux", new ExprFieldWrapper(velocity[0]));
      w.addField("uy", new ExprFieldWrapper(velocity[1]));
      w.addField("uz", new ExprFieldWrapper(velocity[2]));
      w.write();


    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
