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


%rename(NonlinearProblem) makeNonlinearProblem;



%inline %{
  /* create a nonlinear operator */
  TSFExtended::
    NonlinearOperator<double> makeNonlinearProblem(
      const SundanceStdMesh::Mesh& mesh, 
      const SundanceCore::Expr& eqn,
      const SundanceCore::Expr& bc,
      const SundanceCore::Expr& v, 
      const SundanceCore::Expr& u,
      const SundanceCore::Expr& u0,
      const TSFExtended::VectorType<double>& vecType)
  {
    return  new SundanceStdFwk::NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType);
  }

  %}


