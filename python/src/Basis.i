// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceStdFwk
{
  class BasisFamily
  {
  public:
    BasisFamily();
    ~BasisFamily();
  };

  %extend BasisFamily
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

%rename(Lagrange) makeLagrange;
%rename(List) makeList;

%inline %{
  /* Create a line mesher */
  SundanceStdFwk::BasisFamily makeLagrange(int order)
  {
    return new SundanceStdFwk::Lagrange(order);
  }
  %}




