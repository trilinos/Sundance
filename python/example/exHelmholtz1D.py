#! /usr/bin/env python

# @HEADER

# @HEADER

import math
import setpath
import PySundance


from PySundance import *

class LeftPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt) < 1.0e-10);

def main():

  vecType = EpetraVectorType()
  pi = 4.0*math.atan(1.0)
  nx = 40
  mesher  = PartitionedLineMesher(0.0, pi/4.0, nx);
  mesh = mesher.getMesh();
  basis = Lagrange(2)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  dx = Derivative(0);
  
  quad = GaussianQuadrature(2)
  quad6 = GaussianQuadrature(6)
  
  lpp = LeftPointPredicate()
  
  interior = MaximalCellFilter()
  leftPtTest = PositionalCellPredicate(lpp)
  points = DimensionalCellFilter(0)
  leftPt = points.subset(leftPtTest)
  
  bc = EssentialBC(leftPt, v*(u-cos(x)), quad)
  eqn = Integral(interior, (dx*v)*(dx*u) - v*u, quad)
  
  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

  solver = readSolver("../../../tests-std-framework/Problem/bicgstab.xml");

  soln = prob.solve(solver)

  exactSoln = cos(x) + sin(x)

  diff = pow(soln - exactSoln, 2.0)

  print "error = " , math.sqrt(diff.integral(interior, mesh, quad))
  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
