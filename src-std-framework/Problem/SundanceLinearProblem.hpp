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

#ifndef SUNDANCE_LINEARPROBLEM_H
#define SUNDANCE_LINEARPROBLEM_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFVector.hpp"
#include "TSFVectorType.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

  /** 
   * LinearProblem encapsulates all information needed to form
   * a discrete linear problem. 
   */
  class LinearProblem 
    : public TSFExtended::ObjectWithVerbosity<LinearProblem>
  {
  public:
    /** Empty ctor */
    LinearProblem();
    
    /** Construct with a mesh, equation set, bcs, test and unknown funcs,
     * and a vector type. */
    LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
                  const Expr& test, const Expr& unk, 
                  const TSFExtended::VectorType<double>& vecType);

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** */
    LinearProblem(const RefCountPtr<Assembler>& assembler);
#endif

    /** Solve the problem using the specified solver algorithm */
    Expr solve(const LinearSolver<double>& solver) const ;

    /** Return the status of the last solve */
    SolverState<double> solveStatus() const ;

    /** Solve the problem, writing the solution into the given function */
    SolverState<double> solve(const LinearSolver<double>& solver,
                              Expr& soln) const ;

    /** Return the vector on the right-hand side of the linear equation */
    Vector<double> getRHS() const ;

    /** Return the operator on the left-hand side of the equation */
    LinearOperator<double> getOperator() const ;

    /** Return the map from cells and functions to row indices */
    const RefCountPtr<DOFMapBase>& rowMap() const 
    {return assembler_->rowMap();}
    
    /** Return the map from cells and functions to column indices */
    const RefCountPtr<DOFMapBase>& colMap() const 
    {return assembler_->colMap();}
    
    /** Return the set of row indices marked as 
     * essential boundary conditions */
    const RefCountPtr<Set<int> >& bcRows() const 
    {return assembler_->bcRows();}

    /** Form a solution expression out of a vector obtained from a linear
     * solver */
    Expr formSolutionExpr(const Vector<double>& solnVector) const ;


  private:
      
    /** */
    RefCountPtr<Assembler> assembler_;

    /** */
    mutable LinearOperator<double> A_;

    /** */
    mutable Vector<double> rhs_;

    /** */
    mutable RefCountPtr<SolverState<double> > status_;
    
  };
}


#endif
