/*
 * Copyright (2009) Kevin Long
 * Department of Mathematics and Statistics
 * Texas Tech University
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *  
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *                                                                           
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA                                                                       
 * 
 */

#include "Sundance.hpp"

// This test depends on Exodus, so skip it if Expdus hasn't been enabled. 
#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)

using namespace Sundance;

/** 
 * This example program sets up and solves the Laplace 
 * equation \f$-\nabla^2 u=0\f$. See the
 * document GettingStarted.pdf for more information.
 */
int main(int argc, char** argv)
{
  try
  {
    /* command-line options */
    std::string meshFile="plateWithHole3D-1";
    std::string solverFile = "aztec-ml.xml";
    Sundance::setOption("meshFile", meshFile, "mesh file");
    Sundance::setOption("solver", solverFile, 
      "name of XML file for solver");

    /* Initialize */
    Sundance::init(&argc, &argv);

    /* --- Specify vector representation to be used --- */
    VectorType<double> vecType = new EpetraVectorType();

    /* --- Read mesh --- */
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource meshSrc
      =  new ExodusMeshReader(meshFile, meshType);
    Mesh mesh = meshSrc.getMesh();

    /* --- Specification of geometric regions --- */

    /* Region "interior" consists of all maximal-dimension cells */
    CellFilter interior = new MaximalCellFilter();

    /* Identify boundary regions via labels in mesh */
    CellFilter edges = new DimensionalCellFilter(2);

    CellFilter south = edges.labeledSubset(1);
    CellFilter east = edges.labeledSubset(2);
    CellFilter north = edges.labeledSubset(3);
    CellFilter west = edges.labeledSubset(4);
    CellFilter hole = edges.labeledSubset(5);
    CellFilter down = edges.labeledSubset(6);
    CellFilter up = edges.labeledSubset(7);

    /* --- Symbolic equation definition --- */

    /* Test and unknown function */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    /* Gradient operator */
    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr dz = new Derivative(2);
    Expr grad = List(dx, dy, dz);

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad1 = new GaussianQuadrature(1);
    QuadratureFamily quad2 = new GaussianQuadrature(2);


    /** Write the weak form */
    Expr eqn 
      = Integral(interior, (grad*u)*(grad*v), quad1)
      + Integral(east, v, quad1);

    /* Write the essential boundary conditions */
    Expr h = new CellDiameterExpr();
    Expr bc = EssentialBC(west, v*u/h, quad2);

    /* Set up linear problem */
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    /* --- solve the problem --- */

    /* Create the solver as specified by parameters in 
     * an XML file */
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver(solverFile);

    /* Solve! The solution is returned as an Expr containing a 
    * DiscreteFunction */
    Expr soln = prob.solve(solver);

    /* --- Postprocessing --- */

    /* Project the derivative onto the P1 basis */
    DiscreteSpace discSpace(mesh, List(basis, basis, basis), vecType);
    L2Projector proj(discSpace, grad*soln);
    Expr gradU = proj.project();

    /* Write the solution and its projected gradient to a VTK file */
    FieldWriter w = new VTKWriter("LaplaceDemo3D");
    w.addMesh(mesh);
    w.addField("soln", new ExprFieldWrapper(soln[0]));
    w.addField("du_dx", new ExprFieldWrapper(gradU[0]));
    w.addField("du_dy", new ExprFieldWrapper(gradU[1]));
    w.addField("du_dz", new ExprFieldWrapper(gradU[2]));
    w.write();

    /* Check flux balance */
    Expr n = CellNormalExpr(3, "n");
    CellFilter wholeBdry = east+west+north+south+up+down+hole;
    Expr fluxExpr 
      = Integral(wholeBdry, (n*grad)*soln, quad1); 
    double flux = evaluateIntegral(mesh, fluxExpr);
    Out::root() << "numerical flux = " << flux << std::endl;

    /* --- Let's compute a few other quantities, such as the centroid of
     * the mesh:*/

    /* Coordinate functions let us build up functions of position */
    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr z = new CoordExpr(2);

    Expr xCMExpr = Integral(interior, x, quad1);
    Expr yCMExpr = Integral(interior, y, quad1);
    Expr zCMExpr = Integral(interior, z, quad1);
    Expr volExpr = Integral(interior, 1.0, quad1);
    
    double vol = evaluateIntegral(mesh, volExpr);
    double xCM = evaluateIntegral(mesh, xCMExpr)/vol;
    double yCM = evaluateIntegral(mesh, yCMExpr)/vol;
    double zCM = evaluateIntegral(mesh, zCMExpr)/vol;
    Out::root() << "centroid = (" << xCM << ", " << yCM 
              << ", " << zCM << ")" << std::endl;

    /* Next, compute the first Fourier sine coefficient of the solution on the
     * surface of the hole.*/
    Expr r = sqrt(x*x + y*y);
    Expr sinPhi = y/r;

    /* Use a higher-order quadrature rule for these integrals */
    QuadratureFamily quad4 = new GaussianQuadrature(4);

    Expr fourierSin1Expr = Integral(hole, sinPhi*soln, quad4);
    Expr fourierDenomExpr = Integral(hole, sinPhi*sinPhi, quad2);
    double fourierSin1 = evaluateIntegral(mesh, fourierSin1Expr);
    double fourierDenom = evaluateIntegral(mesh, fourierDenomExpr);
    Out::root() << "fourier sin m=1 = " << fourierSin1/fourierDenom << std::endl;

    /* Compute the L2 norm of the solution */
    Expr L2NormExpr = Integral(interior, soln*soln, quad2);     
    double l2Norm_method1 = sqrt(evaluateIntegral(mesh, L2NormExpr));     
    Out::os() << "method #1: ||soln|| = " << l2Norm_method1 << endl;

    /* Use the L2Norm() function to do the same calculation */
    double l2Norm_method2 = L2Norm(mesh, interior, soln, quad2);     
    Out::os() << "method #2: ||soln|| = " << l2Norm_method2 << endl;

    /*
     * Check that the flux is acceptably close to zero. The flux calculation
     * is only O(h) so keep the tolerance loose. This
     * is just a sanity check to ensure the code doesn't get completely 
     * broken after a change to the library. 
     */
    Sundance::passFailTest(fabs(flux), 1.0e-2);
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 

  return Sundance::testStatus();
}

/*
 * All done!
 */

#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy Laplace3D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif

