// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceLinearProblem.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace SundanceStdFwk
{

  class LinearProblem
  {
  public:
    LinearProblem(const SundanceStdMesh::Mesh& mesh, 
                  const SundanceCore::Expr& eqn,
                  const SundanceCore::Expr& bc,
                  const SundanceCore::Expr& v, 
                  const SundanceCore::Expr& u,
                  const TSFExtended::VectorType<double>& vecType);

    TSFExtended::Vector<double> getRHS() const ;

    TSFExtended::LinearOperator<double> getOperator() const ;

    SundanceCore::Expr solve(const TSFExtended::LinearSolver<double>& solver) const ;
  };
}
