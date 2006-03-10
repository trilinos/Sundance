// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceFIATLagrange.hpp"
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
%rename(FIATLagrange) makeFIATLagrange;

%inline %{
  /* Create a Lagrange basis function */
  SundanceStdFwk::BasisFamily makeLagrange(int order)
  {
    return new SundanceStdFwk::Lagrange(order);
  }
  /* Create a Lagrange basis function */
  SundanceStdFwk::BasisFamily makeFIATLagrange(int order)
  {
    #ifdef HAVE_FIAT
    return new SundanceStdFwk::FIATLagrange(order);
    #else
    TEST_FOR_EXCEPTION(true, RuntimeError, "FIATLagrange called, but "
                       "FIAT not enabled. Try reconfiguring with "
                       "--enable-fiat");
    return BasisFamily(); // -Wall
    #endif
  }

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




