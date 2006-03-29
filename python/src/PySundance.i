// -*- c++ -*-

%module PySundance

%feature("autodoc");

%exception 
{
  try
    {
      $action
    }
  catch (std::exception& e)
    {
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
      return NULL;
    }
}


%{
#include "Sundance.hpp"
  %}

%inline %{
  bool passFailTest(double err, double tol)
  {
    return SundanceStdFwk::Sundance::passFailTest(err, tol);
  }
  %}




%include Mesh.i

%include Utils.i

%include ParameterList.i

%include CellFilter.i

%include TSF.i

%include Quadrature.i

%include Basis.i

%include Symbolics.i

%include Integral.i

%include LinearProblem.i

%include NonlinearProblem.i

%include Viz.i

%include Discrete.i

  
