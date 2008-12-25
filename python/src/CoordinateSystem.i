// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceCoordinateSystem.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceStdFwk
{
  

  /* */
  class CoordinateSystem
  {
  public:
    /* */
    Expr jacobian() const ;

    /** */
    bool supportsMeshDimension(int dim) const ;

  };
}

%rename(CartesianCoordinateSystem) makeCartesianCoordinateSystem;
%rename(MeridionalCylindricalCoordinateSystem) makeMeridionalCylindricalCoordinateSystem;
%rename(RadialSphericalCoordinateSystem) makeRadialSphericalCoordinateSystem;

%inline %{
  /* Create a maximal cell filter */
  SundanceStdFwk::CoordinateSystem makeCartesianCoordinateSystem()
  {
    return new SundanceStdFwk::CartesianCoordinateSystem();
  }
  %}
 
%inline %{
  /* Create a maximal cell filter */
  SundanceStdFwk::CoordinateSystem makeMeridionalCylindricalCoordinateSystem()
  {
    return new SundanceStdFwk::MeridionalCylindricalCoordinateSystem();
  }
  %}
 
%inline %{
  /* Create a maximal cell filter */
  SundanceStdFwk::CoordinateSystem makeRadialSphericalCoordinateSystem()
  {
    return new SundanceStdFwk::RadialSphericalCoordinateSystem();
  }
  %}
 
