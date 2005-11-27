// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceExpr.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceParameter.hpp"
#include "SundanceFunctionalEvaluator.hpp"
#include "SundanceStdMathOps.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceCore
{
  class Expr
  {
  public:
    Expr();
    ~Expr();
  };

  %extend Expr
  {
    using namespace std;
    string __str__() 
    {
      string rtn = self->toString(); 
      return rtn;
    }

    double integral(const SundanceStdFwk::CellFilter& domain,
                    const SundanceStdMesh::Mesh& mesh,
                    const SundanceStdFwk::QuadratureFamily& quad)
    {
      Expr I = Integral(domain, *self, quad);
      return SundanceStdFwk::evaluateIntegral(mesh, I);
    }

    Expr __pow__(const double& ex)
    {
      return pow(*self, ex);
    }

    Expr __add__(const SundanceCore::Expr& other) 
    {
      SundanceCore::Expr rtn = self->operator+(other);
      return rtn;
    }

    Expr __sub__(const SundanceCore::Expr& other) 
    {
      SundanceCore::Expr rtn = self->operator-(other);
      return rtn;
    }

    Expr __mul__(const SundanceCore::Expr& other) 
    {
      SundanceCore::Expr rtn = self->operator*(other);
      return rtn;
    }

    Expr __div__(const SundanceCore::Expr& other) 
    {
      SundanceCore::Expr rtn = self->operator/(other);
      return rtn;
    }

    

    /* operations with scalars to the right */

    Expr __add__(const double& other) 
    {
      SundanceCore::Expr rtn = self->operator+(other);
      return rtn;
    }

    Expr __sub__(const double& other) 
    {
      SundanceCore::Expr rtn = self->operator-(other);
      return rtn;
    }

    Expr __mul__(const double& other) 
    {
      SundanceCore::Expr rtn = self->operator*(other);
      return rtn;
    }

    Expr __div__(const double& other) 
    {
      SundanceCore::Expr rtn = self->operator/(other);
      return rtn;
    }

    /* operations with scalars to the left */

    Expr __radd__(const double& other) 
    {
      SundanceCore::Expr rtn = other + *self;
      return rtn;
    }

    Expr __rsub__(const double& other) 
    {
      SundanceCore::Expr rtn = other - *self;
      return rtn;
    }

    Expr __rmul__(const double& other) 
    {
      SundanceCore::Expr rtn = other * (*self);
      return rtn;
    }

    Expr __div__(const double& other) 
    {
      SundanceCore::Expr rtn = other / (*self);
      return rtn;
    }

    /* unary operations */
    Expr __neg__() 
    {
      SundanceCore::Expr rtn = self->operator-();
      return rtn;
    }


    /* list indexing and information */
    
    Expr __getitem__(int i) const
    {
      return self->operator[](i);
    }


    


    
  }

  Expr List(const Expr& a);
  Expr List(const Expr& a, const Expr& b);
  Expr List(const Expr& a, const Expr& b, const Expr& c);

  Expr pow(const Expr& x, const double& a);

  Expr fabs(const Expr& x);
  Expr sign(const Expr& x);

  Expr exp(const Expr& x);
  Expr log(const Expr& x);
  Expr sqrt(const Expr& x);

  Expr sin(const Expr& x);
  Expr cos(const Expr& x);
  Expr tan(const Expr& x);
  Expr asin(const Expr& x);
  Expr acos(const Expr& x);
  Expr atan(const Expr& x);

  Expr sinh(const Expr& x);
  Expr cosh(const Expr& x);
  Expr tanh(const Expr& x);
  Expr asinh(const Expr& x);
  Expr acosh(const Expr& x);
  Expr atanh(const Expr& x);

}



%rename(UnknownFunction) makeUnknownFunction;
%rename(TestFunction) makeTestFunction;
%rename(CoordExpr) makeCoordExpr;
%rename(Derivative) makeDerivative;
%rename(Parameter) makeParameter;

%inline %{
  /* Create an unknown function */
  SundanceCore::Expr makeUnknownFunction(const SundanceStdFwk::BasisFamily& b,
                                         const std::string& name)
  {
    return new SundanceStdFwk::UnknownFunction(b, name);
  }
  %}

%inline %{
  /* Create a test function */
  SundanceCore::Expr makeTestFunction(const SundanceStdFwk::BasisFamily& b,
                                      const std::string& name)
  {
    return new SundanceStdFwk::TestFunction(b, name);
  }
  %}


%inline %{
  /* Create a coordinate expression */
  SundanceCore::Expr makeCoordExpr(int dir)
  {
    return new SundanceCore::CoordExpr(dir);
  }
  %}


%inline %{
  /* Create a differential operator */
  SundanceCore::Expr makeDerivative(int dir)
  {
    return new SundanceCore::Derivative(dir);
  }
  %}


%inline %{
  /* Create a differential operator */
  SundanceCore::Expr makeParameter(const double& val)
  {
    return new SundanceCore::Parameter(val);
  }
  %}


%inline %{
  /* Create a differential operator */
  SundanceCore::Expr makeParameter(const double& val, const std::string& name)
  {
    return new SundanceCore::Parameter(val, name);
  }
  %}

