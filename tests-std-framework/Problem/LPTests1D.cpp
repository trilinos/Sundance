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
#include "SundanceCubicHermite.hpp"
#include "SundanceProblemTesting.hpp"


/** Solves the Poisson equation in 1D with Lagrange(2) basis functions */
class SimplePoisson1DTest : public LP1DTestBase
{
public:
  /** */
  SimplePoisson1DTest() 
    : LP1DTestBase(10) {}

  /** */
  std::string name() const {return "SimplePoisson1D";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);

      return x*(x-2.0);
    }

  /** */
  LinearProblem prob() const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();

      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");
      Expr dx = gradient(1);

      QuadratureFamily quad = new GaussianQuadrature(4);

      Expr eqn = Integral(interior(), (dx*v)*(dx*u), quad)
        + Integral(interior(), 2.0*v, quad);

      Expr bc = EssentialBC(left, v*u, quad);
      
      return LinearProblem(mesh(), eqn, bc, v, u, vecType());
    }
  
};



/** Solves the Poisson equation in 1D with user-specified basis functions */
class Poisson1DTest : public LP1DTestBase
{
public:
  /** */
  Poisson1DTest(int nx, const BasisFamily& basis) 
    : LP1DTestBase(nx), basis_(basis) {}

  /** */
  std::string name() const {return "Poisson1D(nx=" 
      + Teuchos::toString(domain().nx()) + ", order=" 
      + Teuchos::toString(basis_.order()) + ")";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      int p = basis_.order();
      return pow(1.0+x,p);
    }

  /** */
  LinearProblem prob() const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();

      Expr u = new UnknownFunction(basis_, "u");
      Expr v = new TestFunction(basis_, "v");
      Expr dx = gradient(1);
      Expr x = coord(0);

      int p = basis_.order();

      QuadratureFamily quad = new GaussianQuadrature(2*p);

      Expr source;
      if (p>1) source = p*(p-1.0)*pow(1.0+x, p-2.0);
      else source = 0.0;

      Expr ex = exactSoln();
      Expr eqn = Integral(interior(), (dx*v)*(dx*u) + v*source, quad);
      Expr bc = EssentialBC(left, v*(u-ex), quad)
        + EssentialBC(right, v*(u-ex), quad);
      
      
      return LinearProblem(mesh(), eqn, bc, v, u, vecType());
    }

  /** */
  Array<LPTestSpec> specs() const
    {
      /* use a tolerance of 10^-13 times an estimate of 
       * the condition number. Belos-ifpack needs a looser tolerance. */
      int np = MPIComm::world().getNProc();
      int nx = domain().nx();
      double cond =  pow((basis_.order()+1)*np*nx, 2.0);
      double tol = 1.0e-13 * cond;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("belos-ifpack.xml", tol * 100.0),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("belos-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
private:
  BasisFamily basis_;
};



/** Performs L2 projection in 1D with user-specified basis functions */
class Projection1DTest : public LP1DTestBase
{
public:
  /** */
  Projection1DTest(int nx, const BasisFamily& basis) 
    : LP1DTestBase(nx), basis_(basis) {}

  /** */
  std::string name() const {return "Projection1D(nx=" 
      + Teuchos::toString(domain().nx()) + ", order=" 
      + Teuchos::toString(basis_.order()) + ")";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      int p = basis_.order();
      return pow(x,p);
    }

  /** */
  LinearProblem prob() const
    {
      Expr u = new UnknownFunction(basis_, "u");
      Expr v = new TestFunction(basis_, "v");
      Expr x = coord(0);

      int p = basis_.order();

      QuadratureFamily quad = new GaussianQuadrature(2*p);

      Expr ex = exactSoln();
      Expr eqn = Integral(interior(), v*(u-ex), quad);
      Expr bc;
      
      return LinearProblem(mesh(), eqn, bc, v, u, vecType());
    }

  /** */
  Array<LPTestSpec> specs() const
    {
      /* use a tolerance times an estimate of 
       * the condition number. Belos-ifpack needs a looser tolerance. */
      int np = MPIComm::world().getNProc();
      int nx = domain().nx();
      double cond =  10.0*pow(np*nx, 2.0);
      double tol = 1.0e-12 * cond;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("belos-ifpack.xml", tol * 10.0),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
private:
  BasisFamily basis_;
};


/** Solves the Helmholtz equation in 1D with Lagrange(2) basis functions */
class Helmholtz1DTest : public LP1DTestBase
{
public:
  /** */
  Helmholtz1DTest(int nx) 
    : LP1DTestBase(0.0, atan(1.0), nx) {}

  /** */
  std::string name() const {return "Helmholtz1D(nx=" 
      + Teuchos::toString(domain().nx()) + ")";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      return cos(x) + sin(x);
    }

  /** */
  LinearProblem prob() const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();

      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");
      Expr dx = gradient(1);
      Expr x = coord(0);

      QuadratureFamily quad = new GaussianQuadrature(4);

      Expr eqn = Integral(interior(), (dx*v)*(dx*u) - v*u, quad);
      Expr bc = EssentialBC(left, v*(u-cos(x)), quad);
      
      return LinearProblem(mesh(), eqn, bc, v, u, vecType());
    }

  /** */
  Array<LPTestSpec> specs() const
    {
      /* use a tolerance of 10^-10 times an estimate of 
       * the condition number. */
      int np = MPIComm::world().getNProc();
      int nx = domain().nx();
      double cond = pow(3*np*nx, 2.0);
      double tol = 1.0e-8 * cond;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("belos-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
private:
};




/** */
class Coupled1DTest : public LP1DTestBase
{
public:
  /** */
  Coupled1DTest() : LP1DTestBase(128) {}

  /** */
  std::string name() const {return "Coupled1D";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);

      Expr x2 = x*x;
      Expr x3 = x*x2;
      
      Expr uExact = (1.0/120.0)*x2*x3 - 1.0/36.0 * x3 + 7.0/360.0 * x;
      Expr vExact = 1.0/6.0 * x * (x2 - 1.0);

      return List(vExact, uExact);
    }

  /** */
  LinearProblem prob() const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();


      Expr u = new UnknownFunction(new Lagrange(5), "u");
      Expr v = new UnknownFunction(new Lagrange(3), "v");
      Expr du = new TestFunction(new Lagrange(5), "du");
      Expr dv = new TestFunction(new Lagrange(3), "dv");

      Expr dx = gradient(1);
      Expr x = coord(0);


      QuadratureFamily quad = new GaussianQuadrature(10);

      Expr eqn = Integral(interior(), 
        (dx*du)*(dx*u) + du*v + (dx*dv)*(dx*v) + x*dv, 
        quad);

      Expr bc = EssentialBC(left, du*u + dv*v, quad)
        + EssentialBC(right, du*u + dv*v, quad);


      return LinearProblem(mesh(), eqn, bc, 
        List(dv,du), List(v,u), vecType());
      
    }  

  /** */
  Array<LPTestSpec> specs() const
    {
      return tuple(
        LPTestSpec("amesos.xml", 1.0e-10, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", 1.0e-10),
        LPTestSpec("belos-ifpack.xml", 1.0e-10),
        LPTestSpec("bicgstab.xml", 1.0e-9)
        );
    }
};







/** */
class PoissonCubicHermite1DTest : public LP1DTestBase
{
public:
  /** */
  PoissonCubicHermite1DTest() : LP1DTestBase(32) {}

  /** */
  std::string name() const {return "PoissonCubicHermite1D";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      int p = 3;

      return pow(1.0+x,p);
    }

  /** */
  LinearProblem prob() const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();

      Expr u = new UnknownFunction(new CubicHermite(), "u");
      Expr v = new TestFunction(new CubicHermite(), "v");

      Expr dx = gradient(1);
      Expr x = coord(0);

      int p = 3;

      QuadratureFamily quad = new GaussianQuadrature(2*p);
      Expr source = p*(p-1.0)*pow(1.0+x, p-2.0);

      Expr ex = exactSoln();
      Expr eqn = Integral(interior(), (dx*v)*(dx*u) + v*source, quad);
      Expr bc;

      return LinearProblem(mesh(), eqn, bc, v, u, vecType());
    }
  
  /** */
  Array<LPTestSpec> specs() const
    {
      return tuple(
        LPTestSpec("bicgstab.xml", 1.0e-9)
        );
    }
  
};






int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);
    Tabs::showDepth() = false;
    LinearSolveDriver::solveFailureIsFatal() = false;

    LPTestSuite tests;

    tests.registerTest(rcp(new SimplePoisson1DTest()));

    for (int i=2; i<=5; i++)
    {
      tests.registerTest(rcp(new Poisson1DTest(8, new Lagrange(i))));
      tests.registerTest(rcp(new Poisson1DTest(8, new Bernstein(i))));
      tests.registerTest(rcp(new Projection1DTest(8, new Lagrange(i))));
      tests.registerTest(rcp(new Projection1DTest(8, new Bernstein(i))));
    }

    tests.registerTest(rcp(new Coupled1DTest())); 

    tests.registerTest(rcp(new Helmholtz1DTest(32))); 


#ifdef ENABLE_CUBIC_HERMITE
    tests.registerTest(rcp(new PoissonCubicHermite1DTest()));
#endif

    bool pass = tests.run();

    Out::root() << "total test status: " << pass << endl;

    Sundance::passFailTest(pass);
  }
	catch(exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
