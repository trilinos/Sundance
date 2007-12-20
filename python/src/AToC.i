// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceAToCPointLocator.hpp"
#include "SundanceAToCDensitySampler.hpp"
#include "SundanceCToAInterpolator.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

%template(doubleVector) std::vector<double>;
%template(intVector) std::vector<int>;

namespace SundanceStdFwk
{
  

  /* */
  class AToCPointLocator
  {
  public:
    /* */
    AToCPointLocator(const SundanceStdMesh::Mesh& mesh, 
                     const CellFilter& subdomain,
                     const std::vector<int>& nx);

  };

  /* */
  class AToCDensitySampler
  {
  public:
    /** */
    AToCDensitySampler(const AToCPointLocator& locator,
                       const TSFExtended::VectorType<double>& vecType);
    /** */
    AToCDensitySampler(const AToCPointLocator& locator,
                       const std::vector<double>& origin,
                       const std::vector<double>& rotationalAxis,
                       const TSFExtended::VectorType<double>& vecType);

    /** */
    SundanceCore::Expr sample(const std::vector<double>& positions,
                              const double& particleWeight) const ;

    /** */
    SundanceCore::Expr resetCounts() const ;

    /** */
    void addToCounts(const std::vector<double>& positions,
                     const double& particleWeight,
                     SundanceCore::Expr density) const ;
  };


  class CToAInterpolator
  {
  public:
    /** */
    CToAInterpolator(const AToCPointLocator& locator,
                     const SundanceCore::Expr& field);

    /** */
    void interpolate(const std::vector<double>& positions,
      std::vector<double>& results) const ;

    /** */
    void updateField(const SundanceCore::Expr& field) ;
  };
  

  
 
}
