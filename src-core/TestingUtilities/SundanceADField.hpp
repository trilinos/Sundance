/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_ADFIELD_H
#define SUNDANCE_ADFIELD_H

#include "SundanceDefs.hpp"
#include "SundanceADBasis.hpp"

namespace SundanceTesting
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** 
   *
   */
  class ADField
  {
  public:
    /** */
    ADField(){;}

    /** */
    ADField(const ADBasis& basis, const double& coeff);

    /** */
    static Point& evalPoint() {static Point x(0.123, 0.456, 0.789); return x;}


    /** */
    void setCoeff(const double& c);

    /** */
    const ADBasis& basis() const {return basis_;}

    /** */
    ADReal evaluate() const ;

    /** */
    double coeff() const {return *coeff_;}

    double value() const {return evaluate().value();}

    ADReal operator+(const ADReal& x) const ;

    ADReal operator+(const double& x) const ;

    ADReal operator+(const ADField& x) const ;



    ADReal operator-(const ADReal& x) const ;

    ADReal operator-(const double& x) const ;

    ADReal operator-(const ADField& x) const ;



    ADReal operator*(const ADReal& x) const ;

    ADReal operator*(const double& x) const ;

    ADReal operator*(const ADField& x) const ;

    ADReal operator-() const ;

    

    
  private:
    ADBasis basis_;
    RefCountPtr<double> coeff_;
  };


  inline ADReal operator+(const ADReal& x, const ADField& y)
  {
    return y + x;
  }

  inline ADReal operator+(const double& x, const ADField& y)
  {
    return y + x;
  }

  inline ADReal operator-(const ADReal& x, const ADField& y)
  {
    return -y + x;
  }

  inline ADReal operator-(const double& x, const ADField& y)
  {
    return -y + x;
  }

  inline ADReal operator*(const ADReal& x, const ADField& y)
  {
    return y * x;
  }

  inline ADReal operator*(const double& x, const ADField& y)
  {
    return y * x;
  } 
  
  inline ADReal sin(const ADField& x)
  {
    return sin(x.evaluate());
  }
  
  inline ADReal cos(const ADField& x)
  {
    return cos(x.evaluate());
  }
}



#endif
