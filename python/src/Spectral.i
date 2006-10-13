// -*- c++ -*-


%{

#define HAVE_PY_FIAT
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceSpectralBasis.hpp"
#include "SundanceHermiteSpectralBasis.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceCore
{
  class SpectralBasis
  {
  public:
    SpectralBasis();
    ~SpectralBasis();

    int getDim() const ;
    
    int getOrder() const ;

    int nterms() const ;

    int getElement(int i) const ;

    double expectation(int i, int j, int k) const ;

    string toString() const ;
  };

  %extend SpectralBasis
  {
    using namespace std;
    string __str__() 
    {
      return self->toString();
    }
  }


}



%rename(HermiteSpectralBasis) makeHermiteSpectralBasis;
%rename(SpectralExpr) makeSpectralExpr;

%inline %{
  /* Create a Hermite basis */
  SundanceCore::SpectralBasis makeHermiteSpectralBasis(int dim, int order)
  {
    return new SundanceCore::HermiteSpectralBasis(dim, order);
  }

  /* Create a Hermite basis */
  SundanceCore::SpectralBasis makeHermiteSpectralBasis(int dim, int order, int nterms)
  {
    return new SundanceCore::HermiteSpectralBasis(dim, order, nterms);
  }

  /* Create a spectral expression */
  SundanceCore::Expr makeSpectralExpr(const SundanceCore::SpectralBasis& sb, 
                                      const SundanceCore::Expr& coeffs)
  {
    return new SundanceCore::SpectralExpr(sb, coeffs);
  }


  %}




