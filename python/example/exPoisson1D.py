#!/usr/bin/env python

import setpath
import PySundance


import math
from PySundance import *
from aztecSolver import solverParams

class LeftPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt) < 1.0e-10);

tsfSolverParams = ParameterList({
    "Max Iterations" : 1000,
    "Tolerance" : 1.0e-12,
    "Restart" : 1000,
    "Verbosity" : 0
    })

ilukParams = ParameterList({
  "Graph Fill" : 1
    })

def main():
  """Poisson example code"""
  vecType = EpetraVectorType()
  nProc = getNProc()
  mesher  = PartitionedLineMesher(0.0, 1.0, 10*nProc);
  mesh = mesher.getMesh();
  basis = FIATLagrange(2)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  dx = Derivative(0);

  quad = GaussianQuadrature(2)
  
  lpp = LeftPointPredicate()
  
  interior = MaximalCellFilter()
  leftPtTest = PositionalCellPredicate(lpp)
  points = DimensionalCellFilter(0)
  leftPt = points.subset(leftPtTest)
  
  bc = EssentialBC(leftPt, v*u, quad)
  eqn = Integral(interior, -(dx*v)*(dx*u) - 2*v, quad)
  
  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

  precond = ILUKPreconditionerFactory(ilukParams)
  solver = GMRESSolver(tsfSolverParams, precond)

  soln = prob.solve(solver)

  # do a silly vector manipulation just to show how to use vectors
  solnVec = soln.getVector()
  print solnVec
  solnVec = 2.0*solnVec - solnVec
  print solnVec
  soln.setVector(solnVec)

  exactSoln = x*(x-2.0)

  diff = (soln - exactSoln)**2.0
  diffDeriv = (dx*(soln - exactSoln))**2.0

  error = math.sqrt(diff.integral(interior, mesh, quad))
  derivError = math.sqrt(diffDeriv.integral(interior, mesh, quad))
  print "error = " , error
  print "deriv error = " , derivError

  error = max(error, derivError)

  print "max err = ", error

  tol = 1.0e-11
  passFailTest(error, tol)
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
