#include "SundanceCoordinateSystem.hpp"

using namespace SundanceStdFwk;

CoordinateSystem 
CoordinateSystemBuilder::makeCoordinateSystem(const std::string& name)
{
  RefCountPtr<CoordinateSystemBase> rtn;

  if (name=="Cartesian")
  {
    rtn = rcp(new CartesianCoordinateSystem());
  }
  else if (name=="Meridional Cylindrical")
  {
    rtn = rcp(new MeridionalCylindricalCoordinateSystem());
  }
  else if (name=="Radial Spherical")
  {
    rtn = rcp(new RadialSphericalCoordinateSystem());
  }
  else
  {
    TEST_FOR_EXCEPTION(true, RuntimeError,
      "coordinate system type=[" << name << "] not recognized");
  }
  
  return rtn;
}


