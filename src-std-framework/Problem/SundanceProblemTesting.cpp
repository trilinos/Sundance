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

#include "SundanceProblemTesting.hpp"

#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFLinearSolverBuilder.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceCoordExpr.hpp"

namespace Sundance
{

bool checkErrorNorms(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& numSoln,
  const Expr& exactSoln,
  const QuadratureFamily& quad,
  double L2Tol,
  double H1SemiTol,
  double H1Tol)
{
  Tabs tab;
  Tabs tab1;
  Expr df = numSoln - exactSoln;

  double L2Err = L2Norm(mesh, filter, df, quad);
  Out::root() << tab << "L2 norm check:" << endl
              << tab1 << setw(16) << setprecision(6) << L2Err 
              << " tol=" << setw(12) << L2Tol ;
  if (fabs(L2Err) > L2Tol) Out::root() << "   <==== FAIL!" << endl;
  else Out::root() << endl;

  double H1SemiErr = H1Seminorm(mesh, filter, df, quad);
  Out::root() << tab << "H1 seminorm check:" << endl
              << tab1 << setw(16) << setprecision(6) << H1SemiErr 
              << " tol=" << setw(12) << H1SemiTol ;
  if (fabs(H1SemiErr) > H1SemiTol) Out::root() << "   <==== FAIL!" << endl;
  else Out::root() << endl;


  double H1Err = H1Norm(mesh, filter, df, quad);
  Out::root() << tab << "H1 norm check:" << endl
              << tab1 << setw(16) << setprecision(6) << H1Err 
              << " tol=" << setw(12) << H1Tol;
  if (fabs(H1Err) > H1Tol) Out::root() << "   <==== FAIL!" << endl;
  else Out::root() << endl;

  return (fabs(L2Err) <= L2Tol) 
    && (fabs(H1SemiErr) <= H1SemiTol) 
    && (fabs(H1Err) <= H1Tol) ;
}


/* -------- LineDomain ----------------- */

LineDomain::LineDomain(int nx)
  : a_(0.0), b_(1.0), nx_(nx), interior_(new MaximalCellFilter()),
    left_(), right_(), mesh_()
{init();}


LineDomain::LineDomain(double a, double b, int nx)
  : a_(a), b_(b), nx_(nx), interior_(new MaximalCellFilter()),
    left_(), right_(), mesh_()
{init();}

void LineDomain::init()
{
  int np = MPIComm::world().getNProc();
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedLineMesher(a_, b_, nx_*np, meshType);
  mesh_ = mesher.getMesh();
  
  
  CellFilter points = new DimensionalCellFilter(0);
  left_ = points.subset(new CoordinateValueCellPredicate(0,a_));
  right_ = points.subset(new CoordinateValueCellPredicate(0,b_));
}



Array<LPTestSpec> LPTestBase::specs() const
{
  return tuple(
    LPTestSpec("amesos.xml", 1.0e-10, makeSet<int>(1)),
    LPTestSpec("aztec-ifpack.xml", 1.0e-10),
    LPTestSpec("aztec-ml.xml", 1.0e-10),
    LPTestSpec("belos-ifpack.xml", 1.0e-8),
    LPTestSpec("belos-ml.xml", 1.0e-10)
    );
}

/** \relates LPTestSpec */
std::ostream& operator<<(std::ostream& os, const LPTestSpec& spec)
{
  os << "LPTestSpec(tol=" << spec.tol() << ", solver=" << spec.solverFile()
     << endl;
  return os;
}



bool LPTestSuite::run() const 
{
  int np = MPIComm::world().getNProc();

  bool allOK = true;

  for (int i=0; i<tests_.size(); i++)
  {
    Tabs tab(0);
    Array<LPTestSpec> specs = tests_[i]->specs();
    for (int j=0; j<specs.size(); j++)
    {
      if (specs[j].nProcIsAllowed(np))
      {
        Out::root() << endl;
        Out::root() << endl;
        Out::root() << endl;
        Out::root() << tab << "running test " << tests_[i]->name()
                    << " with spec " << specs[j] << endl;

        std::string solverFile = specs[j].solverFile();
        double tol = specs[j].tol();
        bool pass = tests_[i]->run(solverFile, tol);
        allOK = pass && allOK;
      }
      else
      {
        Out::root() << tab << "skipping test " << tests_[i]->name()
                    << " with spec=" << specs[j] << endl;
      }
    }
    Out::root() << endl;
    Out::root() << endl;
  }
  return allOK;
}


LPTestSuite::LPTestSuite()
  : tests_(), testSpecs_() {}


void LPTestSuite::registerTest(const RCP<LPTestBase>& test) 
{
  tests_.append(test);
}

VectorType<double> LPTestBase::vecType() const
{
  return new EpetraVectorType();
}

QuadratureFamily LPTestBase::testQuad() const
{
  return new GaussianQuadrature(8);
}

Expr LPTestBase::coord(int d) const 
{
  TEST_FOR_EXCEPT(d<0 || d>2);
  return new CoordExpr(d);
}

bool LPTestBase::run(const std::string& solverFile, double tol) const
{
  Tabs tab(0);
  LinearProblem lp = prob();

  

  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver(solverFile);

  Expr soln;
  SolverState<double> state = lp.solve(solver, soln);

  if (state.finalState() != SolveConverged) 
  {
    Out::root() << tab << "solver state: " << state << endl;
    Out::root() << tab << "Solve failed!" << endl;
    return false;
  }

  bool allPass = true;
  Expr exact = exactSoln();
  for (int i=0; i<soln.size(); i++)
  {
    Out::root() << tab << "checking solution component i=" << i << endl;
    bool pass = checkErrorNorms(mesh(), interior(), 
      soln[i], exact[i],
      testQuad(), tol, tol, tol);
    allPass = pass && allPass;
  }

  return allPass;
}



LP1DTestBase::LP1DTestBase(int nx)
  : domain_(nx) {}

LP1DTestBase::LP1DTestBase(double a, double b, int nx)
  : domain_(a, b, nx) {}

}
