/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_UNARYFUNCTOR_H
#define SUNDANCE_UNARYFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceExceptions.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     * 
     */
    class UnaryFunctor
    {
    public:
      /** ctor */
      UnaryFunctor(const string& name) : name_(name) {;}

      /** */
      const string& name() const {return name_;}

      /** */
      virtual void eval(const double& x, double* vals, int nDerivs) const = 0 ;
    private:
      string name_;
    };

    /**
     *
     */
    class StdExp : public UnaryFunctor
    {
    public:
      /** */
      StdExp() : UnaryFunctor("exp") {;}
      
      /** */
      void eval(const double& x, double* vals, int nDerivs) const 
      {
        double f = ::exp(x);
        for (int i=0; i<=nDerivs; i++) vals[i] = f;
      }
    };

    /**
     *
     */
    class StdReciprocal : public UnaryFunctor
    {
    public:
      /** */
      StdReciprocal() : UnaryFunctor("reciprocal") {;}
      
      /** */
      void eval(const double& x, double* vals, int nDerivs) const 
      {
        double f = 1.0/x;
        if (nDerivs==0)
          {
            vals[0] = f;
            return;
          }
        if (nDerivs==1)
          {
            vals[0] = f;
            vals[1] = -f*f;
            return;
          }
        TEST_FOR_EXCEPTION(nDerivs > 1, RuntimeError,
                           "StdReciprocal does not implement derivs of order > 1");
        
      }
    };

    /** */
    class StdCos : public UnaryFunctor
    {
    public:
      /** */
      StdCos() : UnaryFunctor("cos") {;}

      /** */
      void eval(const double& x, double* vals, int nDerivs) const 
      {
        double c = ::cos(x);
        if (nDerivs == 0)
          {
            vals[0] = c;
            return;
          }
        double s = ::sin(x);
        if (nDerivs==1) 
          {
            vals[0] = c;
            vals[1] = -s;
            return;
          }
        if (nDerivs==2) 
          {
            vals[0] = c;
            vals[1] = -s;
            vals[2] = -c;
            return;
          }
        if (nDerivs==3) 
          {
            vals[0] = c;
            vals[1] = -s;
            vals[2] = -c;
            vals[3] = s;
            return;
          }
        TEST_FOR_EXCEPTION(nDerivs > 3, RuntimeError,
                           "StdCos does not implement derivs of order > 3");
      }
      
    };
    
    /** */
    class StdSin : public UnaryFunctor
    {
    public:
      /** */
      StdSin() : UnaryFunctor("sin") {;}

      /** */
      void eval(const double& x, double* vals, int nDerivs) const 
      {
        double s = ::sin(x);
        if (nDerivs == 0)
          {
            vals[0] = s;
            return;
          }
        double c = ::cos(x);
        if (nDerivs==1) 
          {
            vals[0] = s;
            vals[1] = c;
            return;
          }
        if (nDerivs==2) 
          {
            vals[0] = s;
            vals[1] = c;
            vals[2] = -s;
            return;
          }
        if (nDerivs==3) 
          {
            vals[0] = s;
            vals[1] = c;
            vals[2] = -s;
            vals[3] = -c;
            return;
          }
        TEST_FOR_EXCEPTION(nDerivs > 3, RuntimeError,
                           "StdSin does not implement derivs of order > 3");
      }
      
    };


  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
