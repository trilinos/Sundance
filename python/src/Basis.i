// -*- c++ -*-


%{


  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#ifdef HAVE_FIAT
#include "SundanceFIATLagrange.hpp"
#include "PySundanceFIATScalarAdapter.hpp"
#include "PySundanceBasisCheck.hpp"
#endif
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


  class BasisArray
  {
  public:
    BasisArray();
    BasisArray(int n);

    void append(const BasisFamily& b);

  };

  %extend BasisArray
  {
    using namespace std;
    string __str__() 
    {
      string rtn; 
      stringstream os;
      os << *self;
      rtn = os.str();
      return rtn;
    }
  }
    

}



%rename(Lagrange) makeLagrange;
#ifdef HAVE_FIAT
%rename(FIATLagrange) makeFIATLagrange;
%rename(FIATScalarAdapter) makeFIATScalarAdapter;
#endif

%inline %{
  /* Create a Lagrange basis function */
  SundanceStdFwk::BasisFamily makeLagrange(int order)
  {
    return new SundanceStdFwk::Lagrange(order);
  }
#ifdef HAVE_FIAT
  SundanceStdFwk::BasisFamily makeFIATLagrange(int order)
  {

    return new SundanceStdFwk::FIATLagrange(order);
  }

  SundanceStdFwk::BasisFamily makeFIATScalarAdapter(PyObject *py_basis ,
						    int order)
  {
    return new SundanceStdFwk::FIATScalarAdapter(py_basis,order);
  }

#endif
  /* */
  SundanceStdFwk::BasisArray 
    BasisList()
  {
    return BasisArray();
  }

  /* */
  SundanceStdFwk::BasisArray 
    BasisList(const SundanceStdFwk::BasisFamily& a)
  {
    return tuple(a);
  }

  /* */
  SundanceStdFwk::BasisArray 
    BasisList(const SundanceStdFwk::BasisFamily& a,
              const SundanceStdFwk::BasisFamily& b)
  {
    return tuple(a,b);
  }

  /* */
  SundanceStdFwk::BasisArray 
    BasisList(const SundanceStdFwk::BasisFamily& a,
              const SundanceStdFwk::BasisFamily& b,
              const SundanceStdFwk::BasisFamily& c)
  {
    return tuple(a,b,c);
  }

  /* */
  SundanceStdFwk::BasisArray 
    BasisList(const SundanceStdFwk::BasisFamily& a,
              const SundanceStdFwk::BasisFamily& b,
              const SundanceStdFwk::BasisFamily& c,
              const SundanceStdFwk::BasisFamily& d)
  {
    return tuple(a,b,c,d);
  }

  /* */
  SundanceStdFwk::BasisArray 
    BasisList(const SundanceStdFwk::BasisFamily& a,
              const SundanceStdFwk::BasisFamily& b,
              const SundanceStdFwk::BasisFamily& c,
              const SundanceStdFwk::BasisFamily& d,
              const SundanceStdFwk::BasisFamily& e)
  {
    return tuple(a,b,c,d,e);
  }
                                      
  %}




