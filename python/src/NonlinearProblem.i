// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceNonlinearProblem.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


 // SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace SundanceStdFwk
{
class NonlinearProblem
{
public:
  /** Empty ctor */
  NonlinearProblem();

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * and a vector type */
  NonlinearProblem(const SundanceStdMesh::Mesh& mesh, const SundanceCore::Expr& eqn, const SundanceCore::Expr& bc,
    const SundanceCore::Expr& test, const SundanceCore::Expr& unk, const SundanceCore::Expr& u0, 
    const TSFExtended::VectorType<double>& vecType,
    bool partitionBCs = false);

  /** Construct with a mesh, equation set, bcs, test and unknown funcs,
   * parameters, and a vector type */
  NonlinearProblem(const SundanceStdMesh::Mesh& mesh, const SundanceCore::Expr& eqn, const SundanceCore::Expr& bc,
    const SundanceCore::Expr& test, const SundanceCore::Expr& unk, const SundanceCore::Expr& u0, 
    const SundanceCore::Expr& params, const SundanceCore::Expr& paramVals,  
    const TSFExtended::VectorType<double>& vecType,
    bool partitionBCs = false);


  

  /** Compute direct sensitivities to parameters */
  SundanceCore::Expr computeSensitivities(const TSFExtended::LinearSolver<double>& solver) const ;

  /** Solve the nonlinear problem */
  NOX::StatusTest::StatusType solve(const TSFExtended::NOXSolver& solver) const ;

  /** Return the current evaluation point as a Sundance expression */
  SundanceCore::Expr getU0() const ;

  /** Set an initial guess */
  void setInitialGuess(const SundanceCore::Expr& u0New);
};

}

