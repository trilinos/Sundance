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


int main(int argc, char** argv)
{
  try
  {
    const double pi = 4.0*atan(1.0);
    double lambda = 1.25*pi*pi;

    int nx = 32;
    int nt = 10;
    double tFinal = 1.0/lambda;

    Sundance::setOption("nx", nx, "Number of elements");
    Sundance::setOption("nt", nt, "Number of timesteps");
    Sundance::setOption("tFinal", tFinal, "Number of timesteps");
    
    Sundance::init(&argc, &argv);

    /* Creation of vector type */
    VectorType<double> vecType = new EpetraVectorType();

    /* Set up mesh */
    MeshType meshType = new BasicSimplicialMeshType();
      
    MeshSource meshSrc = new PartitionedRectangleMesher(
      0.0, 1.0, nx,
      0.0, 1.0, nx,
      meshType);
    Mesh mesh = meshSrc.getMesh();

    /* 
     * Specification of cell filters
     */
    CellFilter interior = new MaximalCellFilter();
    CellFilter edges = new DimensionalCellFilter(1);
    CellFilter west = edges.coordSubset(0, 0.0);
    CellFilter east = edges.coordSubset(0, 1.0);
    CellFilter south = edges.coordSubset(1, 0.0);
    CellFilter north = edges.coordSubset(1, 1.0);

    /* set up test and unknown functions */
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    /* set up differential operators */
    Expr grad = gradient(2);

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);


    DiscreteSpace discSpace(mesh, basis, vecType);
    Expr uExact = cos(0.5*pi*y)*sin(pi*x)*exp(-lambda*t);
    L2Projector proj(discSpace, uExact);
    Expr uPrev = proj.project();


    /* 
     * We need a quadrature rule for doing the integrations 
     */
    QuadratureFamily quad = new GaussianQuadrature(2);

    double deltaT = tFinal/nt;

    Expr gWest = -pi*exp(-lambda*t)*cos(0.5*pi*y);
    Expr gWestPrev = -pi*exp(-lambda*tPrev)*cos(0.5*pi*y);
    
    /* Create the weak form */
    Expr eqn = Integral(interior, v*(u-uPrev)/deltaT
      + 0.5*(grad*v)*(grad*u + grad*uPrev), quad)
      + Integral(west, -0.5*v*(gWest+gWestPrev), quad);

    Expr bc = EssentialBC(east + north, v*u, quad);

    
    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    
    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    FieldWriter w0 = new VTKWriter("TransientHeat2D-0");
    w0.addMesh(mesh);
    w0.addField("T", new ExprFieldWrapper(uPrev[0]));
    w0.write();

    for (int i=0; i<nt; i++)
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


    
    double err = L2Norm(mesh, interior, uExact-uPrev, quad);
    Out::root() << "error norm=" << err << endl;

    double h = 1.0/(nx-1.0);
    double tol = 0.1*(pow(h,2.0) + pow(lambda*deltaT, 2.0));
    Out::root() << "tol=" << tol << endl;
    
    
    Sundance::passFailTest(err, tol);
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus();
}

