#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

######################################################################
#
# Test of a nonlinear problem where some functions are defined on
# subdomains only
#
######################################################################


class APredicate :
  def evalOp(self, x, y) :
    return (x <= 0.4);

class BPredicate :
  def evalOp(self, x, y) :
    return (x >= 0.4 and x <= 0.6);

from noxSolver import solverParams

def main():

  vecType = EpetraVectorType()
  nx = 200
  ny = 10
  mesher  = PartitionedRectangleMesher(0.0, 1.0, nx, 1,
                                       0.0, 2.0, ny, 1);
  mesh = mesher.getMesh();
  basis = Lagrange(1)

  
  u1 = UnknownFunction(basis);
  v1 = TestFunction(basis);
  u2 = UnknownFunction(basis);
  v2 = TestFunction(basis);
  
  x = CoordExpr(0);
  
  quad = GaussianQuadrature(4)

  # define cell filters for the entire domain and the subsets
  # (0.0, 0.4) and (0.4, 0.6)
  interior = MaximalCellFilter()
  A = interior.subset(APredicate())
  B = interior.subset(BPredicate())

  # Now we define a discrete space for two field variables, one of
  # which is defined only on a subset of the domain.
  # To specify the domains of the functions, create a CellFilterList
  # containing the domains for each function. In the present
  # example, the first field is discretized on the whole interior
  # and the second on the union of A and B. 
  funcDomains = CellFilterList(interior, A+B)
  basisList = BasisList(basis, basis)
  discSpace = DiscreteSpace(mesh, basisList, funcDomains, vecType)
  # Define a discrete function on this space. Note that projections
  # onto partial-domain discrete spaces are not yet available. 
  u0 = DiscreteFunction(discSpace, 1.0)

  # set up an equation
  eqn = Integral(interior, v1*(u1 - 2.0), quad) \
        + Integral(A, v2*(u2*u2 - x*u1), quad) \
        + Integral(B, v2*(u2*u2 - 0.4*u1), quad)

  bc = Expr()

  # set up the problem
  prob = NonlinearProblem(mesh, eqn, bc,
                          List(v1, v2), List(u1, u2), u0, vecType)

  # solve the problem
  solver = NOXSolver(solverParams, prob)
  solver.solve()

  # compute the error
  err2 = ((u0[0] - 2.0)**2.0).integral(interior, mesh, quad) \
        + ((u0[1]*u0[1] - x*u0[0])**2.0).integral(A, mesh, quad)\
        + ((u0[1]*u0[1] - 0.4*u0[0])**2.0).integral(B, mesh, quad)

  error = math.sqrt(err2)
  print "error = " , error

  writer = VTKWriter("PartialDomain2D");
  writer.addMesh(mesh)
  writer.addField("u0", u0[0])
  writer.addField("u1", u0[1])
  writer.write()

  tol = 1.0e-4
  passFailTest(error, tol)
  
  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
