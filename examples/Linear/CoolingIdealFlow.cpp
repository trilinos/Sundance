/* 
 * <Ignore> 
 * Copyright (2012) Kevin Long
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
using Sundance::List;

// This test depends on Exodus, so skip it if Expdus hasn't been enabled. 
#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)


/* 
 * This example simulates heat transfer in an ideal fluid passing by 
 * a solid cylinder. Transport in the fluid is advection and diffusion.
 * Transport in the solid cylinder is diffusion only.
 */


int main(int argc, char** argv)
{
  try
  {
    std::string meshFile="pipe2D-2";
    std::string outFile="TwoZoneTransport1";
    double K_f = 0.145; // conductivity of fluid (W/cm/K)
    double c_f = 2.0; // specific heat of fluid  (J/cm/K)
    double rho_f = 0.865; // density fluid  (g/cm^3)
    double K_s = 0.61; // conductivity of solid (W/cm/K)
    double U0 = 20.0; // average velocity (cm/s)
    double Q = 1.0; // heat generation rate in solid (W/cm^3)
    std::string heatSolverFile="aztec-ifpack.xml";
    std::string flowSolverFile="aztec-ml.xml";

    Sundance::setOption("meshFile", meshFile, "mesh file (omit .exo suffix)");
    Sundance::setOption("outFile", outFile, "output file (omit .vtu suffix)");
    Sundance::setOption("K-fluid", K_f, "Thermal diffusivity of fluid");
    Sundance::setOption("rho-fluid", rho_f, "Density of fluid");
    Sundance::setOption("c-fluid", rho_f, "Specific heat of fluid");
    Sundance::setOption("K-solid", K_s, "Thermal diffusivity of solid");
    Sundance::setOption("U0", U0, "Average velocity of fluid");
    Sundance::setOption("Q", Q, "Heat generation rate in solid");
    Sundance::setOption("flow-solver", flowSolverFile, "flow solver XML file");
    Sundance::setOption("heat-solver", heatSolverFile, "heat solver XML file");

    Sundance::init(&argc, &argv);

    Out::root() << "K_f = " << K_f << endl;
    Out::root() << "rho_f = " << rho_f << endl;
    Out::root() << "c_f = " << c_f << endl;
    Out::root() << "K_s = " << K_s << endl;
    Out::root() << "Q = " << Q << endl;
    Out::root() << "Peclet number = " << rho_f*c_f*U0 / K_f << endl;

    /* use epetra */
    VectorType<double> vecType = new EpetraVectorType();

    /* read solver parameters */
    LinearSolver<double> flowSolver 
      = LinearSolverBuilder::createSolver(flowSolverFile);

    LinearSolver<double> heatSolver 
      = LinearSolverBuilder::createSolver(heatSolverFile);


    /* get mesh */
    Out::root() << "reading mesh" << endl;      
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource meshSrc
      =  new ExodusMeshReader(meshFile, meshType);
    Mesh mesh = meshSrc.getMesh();

    /* set up cell filters */
    CellFilter interior = new MaximalCellFilter();
    CellFilter fluid = interior.labeledSubset(1);
    CellFilter solid = interior.labeledSubset(2);

    CellFilter sides = new DimensionalCellFilter(1);

    CellFilter east = sides.labeledSubset(1);
    CellFilter west = sides.labeledSubset(2);
    CellFilter north = sides.labeledSubset(3);
    CellFilter south = sides.labeledSubset(4);


    /* */
    BasisFamily basis = new Lagrange(1);
    Expr phi = new UnknownFunction(basis, "phi");
    Expr v = new TestFunction(basis, "v");

    /* */
    Expr grad = gradient(2);

    /* */
    QuadratureFamily quad2 = new GaussianQuadrature(2);
    QuadratureFamily quad4 = new GaussianQuadrature(4);

    /* Potential flow equations */
    Expr flowEqn = Integral(fluid, (grad*phi)*(grad*v), quad2);

    Expr h = new CellDiameterExpr();
    Expr x = new CoordExpr(0);
    Expr flowBC = EssentialBC(west, v*(phi-x*U0)/h, quad2) +
      EssentialBC(east, v*(phi-x*U0)/h, quad2);

    /* Flow problem */
    LinearProblem flowProb(mesh, flowEqn, flowBC, v, phi, vecType);

    Out::root() << "solving flow" << endl;
    Expr phi0 = flowProb.solve(flowSolver);

    /* Heat transport equations */
    Expr T = new UnknownFunction(basis);

    Expr n = CellNormalExpr(2, "n");

    Expr heatEqn = Integral(solid, K_s*(grad*v)*(grad*T) - v*Q, quad2)
      + Integral(fluid, 
        (grad*v)*(K_f*(grad*T) - c_f*rho_f*(grad*phi0)*T), quad2)
      + Integral(east, c_f*rho_f*v*T*n*(grad*phi0), quad2);

    Expr heatBC = EssentialBC(west, v*T/h, quad2);

    /* Flow problem */
    LinearProblem heatProb(mesh, heatEqn, heatBC, v, T, vecType);

    Out::root() << "solving heat" << endl;
    Expr T0 = heatProb.solve(heatSolver);

    Out::root() << "writing output" << endl;
    FieldWriter w = new VTKWriter(outFile);
    w.addMesh(mesh);
    w.addField("phi", new ExprFieldWrapper(phi0[0]));
    w.addField("T", new ExprFieldWrapper(T0[0]));
    w.write();


    /* Validation: check that energy is conserved */
    Out::root() << "checking flux balance" << endl;
    Expr J = c_f*rho_f*T0*(grad*phi0) - K_f*(grad*T0);
    Expr fluxCheckNet = Integral(east+west+north+south, n*J, quad2);
    Expr sourceCheck = Integral(solid, Q, quad2);

    double netFlux = evaluateIntegral(mesh, fluxCheckNet);
    double netSource = evaluateIntegral(mesh, sourceCheck);

    
    Out::root() << endl;
    Out::root() << "net heat flux out: " << netFlux << endl;
    Out::root() << "net heat production by source: " << netSource << endl;
    Out::root() << endl;


    

    Sundance::passFailTest(fabs(netFlux - netSource), 1.0e-2);
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 

  return Sundance::testStatus();
}

#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy PoissonDemo3D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif


