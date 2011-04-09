/* 
 * <Ignore> 
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
 *</Ignore>
 */

#include "Sundance.hpp"


CELL_PREDICATE(Omega1Test, {return x[0] <= 0.5;})
CELL_PREDICATE(Omega2Test, {return x[0] >= 0.5;})

int main(int argc, char** argv)
{
  try
  {
    /*
     * Initialization code
     */
    Sundance::init(&argc, &argv);

    /* Creation of vector type */
    VectorType<double> vecType = new EpetraVectorType();

    /* Set up mesh */
    int nx = 32;
    int npx = 1;
    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource meshSrc = new PartitionedRectangleMesher(
      0.0, 1.0, nx, npx,
      0.0, 1.0, nx, npx,
      meshType);
    Mesh mesh = meshSrc.getMesh();

    /* 
     * Specification of cell filters
     */
    CellFilter interior = new MaximalCellFilter();
    CellFilter omega1 = interior.subset(new Omega1Test());
    CellFilter omega2 = interior.subset(new Omega2Test());
    CellFilter edges = new DimensionalCellFilter(1);
    CellFilter south = edges.subset(new CoordinateValueCellPredicate(1,0.0));
    CellFilter north = edges.subset(new CoordinateValueCellPredicate(1,1.0));

    /* set up test and unknown functions */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    /* set up differential operators */
    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr grad = List(dx, dy);

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    DiscreteSpace discSpace(mesh, basis, vecType);
    Expr uPrev = new DiscreteFunction(discSpace, 300.0, "uPrev");

    double C1 = 0.897;
    double C2 = 0.385;
    double rho1 = 2.7;
    double rho2 = 8.96;
    double kappa1 = 0.237;
    double kappa2 = 0.401;

    const double pi = 4.0*atan(1.0);

    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);

    /* 
     * We need a quadrature rule for doing the integrations 
     */
    QuadratureFamily quad = new GaussianQuadrature(2);

    int nSteps = 32;
    double deltaT = 5.0/nSteps;

    
    /* Create the weak form */
    Expr eqn = Integral(omega1, rho1*C1*v*(u-uPrev)/deltaT
      + 0.5*kappa1*(grad*v)*(grad*u + grad*uPrev), quad)
      + Integral(omega2, rho2*C2*v*(u-uPrev)/deltaT
        + 0.5*kappa2*(grad*v)*(grad*u + grad*uPrev), quad)
      + Integral(north, v*(u-(300.0+20.0*sin(pi*t/2.0))), 
        quad);

    Expr bc = EssentialBC(south, v*(u-300.0), quad);

    
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    for (int i=0; i<nSteps; i++)
    {
      t.setParameterValue((i+1)*deltaT);
      tPrev.setParameterValue(i*deltaT);
      Out::root() << "t=" << (i+1)*deltaT << endl;
      Expr uNext = prob.solve(solver);
      
      ostringstream oss;
      oss << "TransientHeat2D-" << i+1;
      FieldWriter w = new VTKWriter(oss.str());
      w.addMesh(mesh);
      w.addField("T", new ExprFieldWrapper(uNext[0]));
      w.write();

      updateDiscreteFunction(uNext, uPrev);
    }

    
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
    return -1;
  }
  Sundance::finalize(); 

  return 0;
}

