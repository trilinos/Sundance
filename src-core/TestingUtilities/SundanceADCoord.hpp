/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ADCOORD_H
#define SUNDANCE_ADCOORD_H

#include "SundanceDefs.hpp"
#include "SundanceADField.hpp"

namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   *
   */
  class ADCoord
  {
  public:
    /** */
    ADCoord(int dir);

    /** */
    ADReal evaluate() const ;

    double value() const {return evaluate().value();}

    ADReal operator+(const ADReal& x) const ;

    ADReal operator+(const double& x) const ;

    ADReal operator+(const ADCoord& x) const ;



    ADReal operator-(const ADReal& x) const ;

    ADReal operator-(const double& x) const ;

    ADReal operator-(const ADCoord& x) const ;



    ADReal operator*(const ADReal& x) const ;

    ADReal operator*(const double& x) const ;

    ADReal operator*(const ADCoord& x) const ;


    ADReal operator+(const ADField& x) const ;

    ADReal operator-(const ADField& x) const ;

    ADReal operator*(const ADField& x) const ;


    ADReal operator-() const ;

    

    
  private:
    int dir_;
  };


  inline ADReal operator+(const ADReal& x, const ADCoord& y)
  {
    return y + x;
  }

  inline ADReal operator+(const ADField& x, const ADCoord& y)
  {
    return y + x;
  }

  inline ADReal operator+(const double& x, const ADCoord& y)
  {
    return y + x;
  }





  inline ADReal operator-(const ADReal& x, const ADCoord& y)
  {
    return -y + x;
  }

  inline ADReal operator-(const double& x, const ADCoord& y)
  {
    return -y + x;
  }

  inline ADReal operator-(const ADField& x, const ADCoord& y)
  {
    return -y + x;
  }





  inline ADReal operator*(const ADReal& x, const ADCoord& y)
  {
    return y * x;
  }

  inline ADReal operator*(const double& x, const ADCoord& y)
  {
    return y * x;
  } 
  

  inline ADReal operator*(const ADField& x, const ADCoord& y)
  {
    return y * x;
  } 
  

  inline ADReal sin(const ADCoord& x)
  {
    return sin(x.evaluate());
  }

  inline ADReal cos(const ADCoord& x)
  {
    return cos(x.evaluate());
  }

}



#endif
