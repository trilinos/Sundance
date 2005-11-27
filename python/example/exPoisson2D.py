#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

class LeftPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(x) < 1.0e-10);

class RightPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(x-1.0) < 1.0e-10);

class BottomPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(y) < 1.0e-10);

class TopPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(y-2.0) < 1.0e-10);

def main():

  vecType = EpetraVectorType()
  npx = 1
  npy = 1
  n = 10
  mesher  = PartitionedRectangleMesher(0.0, 1.0, n, npx,
                                       0.0, 2.0, n, npy);
  mesh = mesher.getMesh();
  basis = Lagrange(2)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  y = CoordExpr(1);
  dx = Derivative(0);
  dy = Derivative(1);
  grad = List(dx, dy);
  
  quad2 = GaussianQuadrature(2)
  quad4 = GaussianQuadrature(4)
  
  interior = MaximalCellFilter()
  edges = DimensionalCellFilter(1)
  left = edges.subset(PositionalCellPredicate(LeftPointPredicate()))
  right = edges.subset(PositionalCellPredicate(RightPointPredicate()))
  top = edges.subset(PositionalCellPredicate(TopPointPredicate()))
  bottom = edges.subset(PositionalCellPredicate(BottomPointPredicate()))

  one = 1.0
  eqn = Integral(interior, (grad*v)*(grad*u) + one*v, quad2)\
        + Integral(top, -v/3.0, quad2)\
        + Integral(right, -v*(1.5 + (1.0/3.0)*y - u), quad4)
  
  bc = EssentialBC(bottom, v*(u-0.5*x*x), quad4)
  
  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

  solver = readSolver("../../../tests-std-framework/Problem/aztec.xml");

  soln = prob.solve(solver)

  exactSoln = 0.5*x*x + (1.0/3.0)*y;

  diff = (soln - exactSoln)**2.0
  diffDeriv = (grad*(soln - exactSoln))**2.0

  print "error = " , math.sqrt(diff.integral(interior, mesh, quad4))
  print "deriv error = " , math.sqrt(diffDeriv.integral(interior, mesh, quad4))


  writer = VTKWriter("Poisson2D");
  writer.addMesh(mesh)
  writer.addField("u0", soln)
  writer.write()





  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
