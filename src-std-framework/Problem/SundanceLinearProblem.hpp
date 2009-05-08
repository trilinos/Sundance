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
#include "SundanceBlock.hpp"
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
  using namespace SundanceCore;
  using namespace Teuchos;

namespace Internal
{
class Assembler;
}

  /** 
   * LinearProblem encapsulates all information needed to form
   * a discrete linear problem. 
   */
  class LinearProblem 
    : public TSFExtended::ParameterControlledObjectWithVerbosity<LinearProblem>
  {
  public:
    /** Empty ctor */
    LinearProblem();
    
    /** Construct with a mesh, equation set, bcs, test and unknown funcs,
     * and a vector type. */
    LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
      const Expr& test, const Expr& unk, 
      const TSFExtended::VectorType<double>& vecType,
      const ParameterList& verbParams = *defaultVerbParams(),
      bool partitionBCs = false
      );
    
    /** Construct with a mesh, equation set, bcs, and blocks of variables */
    LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
      const BlockArray& test, const BlockArray& unk,
      const ParameterList& verbParams = *defaultVerbParams(),
      bool partitionBCs = false);
    
    /** Construct with a mesh, equation set, bcs, test and unknown funcs,
     * parameters, and a vector type. */
    LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
      const Expr& test, const Expr& unk, 
      const Expr& unkParams, const Expr& unkParamVals, 
      const TSFExtended::VectorType<double>& vecType,
      const ParameterList& verbParams = *defaultVerbParams(),
      bool partitionBCs = false);
    
    /** Construct with a mesh, equation set, bcs, parameters, and blocks of variables */
    LinearProblem(const Mesh& mesh, const Expr& eqn, const Expr& bc,
      const BlockArray& test, const BlockArray& unk, 
      const Expr& unkParams, const Expr& unkParamVals,
      const ParameterList& verbParams = *defaultVerbParams(),
      bool partitionBCs = false);

    /** */
    LinearProblem(const RefCountPtr<Assembler>& assembler,
      const ParameterList& verbParams = *defaultVerbParams());

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
    const RefCountPtr<DOFMapBase>& rowMap(int blockRow) const ;
    
    /** Return the map from cells and functions to column indices */
    const RefCountPtr<DOFMapBase>& colMap(int blockCol) const ;

    /** Return the discrete space in which solutions live */
    const Array<RefCountPtr<DiscreteSpace> >& solnSpace() const ;

    
    /** Return the set of row indices marked as 
     * essential boundary conditions */
    const RefCountPtr<Set<int> >& bcRows(int blockRow) const ;

    /** Return the number of block rows in the problem  */
    int numBlockRows() const ;

    /** Return the number of block cols in the problem  */
    int numBlockCols() const ;

    /** Form a solution expression out of a vector obtained from a linear
     * solver */
    Expr formSolutionExpr(const Vector<double>& solnVector) const ;

    /** Convert from a BC-partitioned solution vector to a 
     * monolithic vector */
    Vector<double> 
    convertToMonolithicVector(const Array<Vector<double> >& internalBlock,
      const Array<Vector<double> >& bcBlock) const ;

    /** Flag indicating whether to stop on a solve failure */
    static bool& stopOnSolveFailure() {static bool rtn = false; return rtn;}

    /** Flag indicating whether to write out the matrix and vector
     * after a solve failure */
    static bool& dumpBadMatrix() {static bool rtn = true; return rtn;}

    /** Filename for dump of bad matrix */
    static string& badMatrixFilename() 
    {static string rtn = "badMatrix.dat"; return rtn;}

    /** Filename for dump of bad vector */
    static string& badVectorFilename() 
    {static string rtn = "badVector.dat"; return rtn;}

    
    /** */
    static RefCountPtr<ParameterList> defaultVerbParams();


  private:

    /** Do error reporting and matrix dumping after a solve failure */
    void handleSolveFailure() const ;
    
      
    /** */
    RefCountPtr<Assembler> assembler_;

    /** */
    mutable LinearOperator<double> A_;

    /** */
    mutable Vector<double> rhs_;

    /** */
    mutable RefCountPtr<SolverState<double> > status_;

    /** */
    Array<Array<string> > names_;
    
  };
}


#endif
