// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceIntegral.hpp"
#include "SundanceEssentialBC.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

%rename(Integral) makeIntegral;
%rename(EssentialBC) makeEssentialBC;


%inline %{
  /* */
  SundanceCore::Expr makeIntegral(const SundanceStdFwk::CellFilter& domain,
                                  const SundanceCore::Expr& integrand,
                                  const SundanceStdFwk::QuadratureFamily& quad)
  {
    return SundanceCore::Integral(domain, integrand, quad);
  }
  %}



%inline %{
  /* */
  SundanceCore::Expr makeEssentialBC(const SundanceStdFwk::CellFilter& domain,
                                     const SundanceCore::Expr& integrand,
                                     const SundanceStdFwk::QuadratureFamily& quad)
  {
    return SundanceCore::EssentialBC(domain, integrand, quad);
  }
  %}


