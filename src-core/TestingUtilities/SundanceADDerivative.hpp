/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ADDERIVATIVE_H
#define SUNDANCE_ADDERIVATIVE_H

#include "SundanceDefs.hpp"
#include "SundanceADField.hpp"
#include "SundanceADCoord.hpp"


namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   *
   */
  class ADDerivative
  {
  public:
    /** */
    ADDerivative(int dir);

    /** */
    ADReal operator*(const double& f) const 
    {return ADReal(0.0, Point(0.0, 0.0, 0.0));}

    /** */
    ADReal operator*(const ADReal& f) const 
    {return ADReal(f.gradient()[dir_], Point(0.0, 0.0, 0.0));}

    /** */
    ADReal operator*(const ADCoord& f) const 
    {return ADReal(f.evaluate().gradient()[dir_], Point(0.0, 0.0, 0.0));}

    /** */
    ADReal operator*(const ADField& f) const  
    {return ADReal(f.evaluate().gradient()[dir_], Point(0.0, 0.0, 0.0));}
    
  private:
    int dir_;
  };



}



#endif
