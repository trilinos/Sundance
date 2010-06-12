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

#ifndef SUNDANCE_PROBLEMTESTING_H
#define SUNDANCE_PROBLEMTESTING_H

#include "SundanceFunctional.hpp"
#include "SundanceLinearProblem.hpp"

namespace Sundance
{

using namespace Teuchos;


/** 
 * This function checks the L2 and H1 norms and the H1 seminorm of
 * an error against a specified tolerance.
 */
bool checkErrorNorms(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& numSoln,
  const Expr& exactSoln,
  const QuadratureFamily& quad,
  double L2Tol,
  double H1SemiTol,
  double H1Tol);



/** 
 * This class bundles together a mesh of the 1D interval [a,b] with cell
 * filters defining the interior and boundaries. It is intended for quick 
 * and reliable setup of 1D test problems.
 */
class LineDomain
{
public:
  /** */
  LineDomain(int nx);

  /** */
  LineDomain(double a, double b, int nx);

  /** */
  const CellFilter& left() const {return left_;}

  /** */
  const CellFilter& right() const {return right_;}

  /** */
  const CellFilter& interior() const {return interior_;}

  /** */
  const Mesh& mesh() const {return mesh_;}

  /** */
  double a() const {return a_;}

  /** */
  double b() const {return b_;}
  
  /** */
  int nx() const {return nx_;}

private:
  void init();

  double a_;
  double b_;
  int nx_;
  CellFilter interior_;
  CellFilter left_;
  CellFilter right_;
  Mesh mesh_;
};



/** */
class LPTestSpec
{
public:
  /** */
  LPTestSpec() {;}
  /** */
  LPTestSpec(const std::string& solverFile, double tol)
    : hasProcRestriction_(false), allowedProcNumbers_(),
      solverFile_(solverFile), tol_(tol){}

  /** */
  LPTestSpec(const std::string& solverFile, double tol, 
    const Set<int>& allowedProcs) 
    : hasProcRestriction_(true), allowedProcNumbers_(allowedProcs),
      solverFile_(solverFile), tol_(tol){}

  /** */
  const double& tol() const {return tol_;}

  /** */
  const std::string& solverFile() const {return solverFile_;}

  /** */
  bool nProcIsAllowed(int np) const
    {
      if (!hasProcRestriction_) return true;
      return allowedProcNumbers_.contains(np);
    }

private:
  bool hasProcRestriction_;

  Set<int> allowedProcNumbers_;

  std::string solverFile_;

  double tol_;
};

/** \relates LPTestSpec */
std::ostream& operator<<(std::ostream& os, const LPTestSpec& spec);

/** */
class LPTestBase
{
public:
  /** */
  virtual bool run(const std::string& solverFile, double tol) const ;

  /** */
  virtual string name() const = 0 ;

  /** */
  virtual Array<LPTestSpec> specs() const ;

  /** */
  virtual Expr exactSoln() const = 0 ;

  /** */
  virtual VectorType<double> vecType() const ;

  /** */
  virtual LinearProblem prob() const = 0 ;

  /** */
  virtual QuadratureFamily testQuad() const ;

  /** */
  virtual Expr coord(int d) const ;

  /** */
  virtual Mesh mesh() const = 0 ;

  /** */
  virtual CellFilter interior() const = 0 ;
};



/** */
class LP1DTestBase : public LPTestBase
{
public:
  /** */
  LP1DTestBase(int nx);

  /** */
  LP1DTestBase(double a, double b, int nx);

  /** */
  CellFilter interior() const {return domain_.interior();}

  /** */
  Mesh mesh() const {return domain_.mesh();}

  /** */
  const LineDomain& domain() const {return domain_;}
  
private:
  LineDomain domain_;
};


/** */
class LPTestSuite
{
public:
  /** */
  LPTestSuite();

  /** */
  void registerTest(const RCP<LPTestBase>& test) ;

  /** */
  bool run() const ;
  
private:
  Array<RCP<LPTestBase> > tests_;
  Array<Array<LPTestSpec> > testSpecs_;
};


}

#endif

