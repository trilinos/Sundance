// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceQuadratureFamily.hpp"
#include "SundanceGaussianQuadrature.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceStdFwk
{
  class QuadratureFamily
  {
  public:
    QuadratureFamily();
    ~QuadratureFamily();
  };

  %extend QuadratureFamily
  {
    using namespace std;
    string __str__() 
    {
      string rtn; 
      stringstream os;
      self->print(os);
      rtn = os.str();
      return rtn;
    }
  }
}

%rename(GaussianQuadrature) makeGaussianQuadrature;
%rename(FIATQuadratureAdapter) makeFIATQuadratureAdapter;

%inline %{
  /* Create a Gaussian quadrature object */
  SundanceStdFwk::QuadratureFamily makeGaussianQuadrature(int order)
  {
    return new SundanceStdFwk::GaussianQuadrature(order);
  }

  /*  SundanceStdFwk::FIATQuadratureAdapter 
  	makeFIATQuadratureAdapter(PyObject *py_quad_factory , int order)
  {
    return new SundanceStdFwk::FIATQuadratureAdapter(factor,order);
    }*/
  
  %}

