/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCTORDOMAIN_H
#define SUNDANCE_FUNCTORDOMAIN_H

#include "SundanceDefs.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  using namespace Teuchos;

  class FunctorDomain
  {
  public:
    FunctorDomain();

    virtual bool hasLowerBound() const {return false;}

    virtual double lowerBound() const ;

    virtual bool hasUpperBound() const {return false;}

    virtual double upperBound() const ;

    virtual bool hasExcludedPoint() const {return false;}

    virtual double excludedPoint() const ;

  };

  class UnboundedDomain : public FunctorDomain
  {
  public:
    UnboundedDomain();
  };


  class PositiveDomain : public FunctorDomain
  {
  public:
    PositiveDomain();

    virtual bool hasLowerBound() const {return true;}

    virtual double lowerBound() const {return 0.0;}
  };


  class BoundedDomain : public FunctorDomain
  {
  public:
    BoundedDomain(const double& lower, const double& upper);

    virtual bool hasLowerBound() const {return true;}

    virtual double lowerBound() const {return lower_;}

    virtual bool hasUpperBound() const {return true;}

    virtual double upperBound() const {return upper_;}

  private:
    double lower_;

    double upper_;
  };


  class LowerBoundedDomain : public FunctorDomain
  {
  public:
    LowerBoundedDomain(const double& lower);

    virtual bool hasLowerBound() const {return true;}

    virtual double lowerBound() const {return lower_;}

  private:
    double lower_;
  };

class NonzeroDomain : public FunctorDomain
  {
  public:
    NonzeroDomain();

    virtual bool hasExcludedPoint() const {return true;}

    virtual double excludedPoint() const {return 0.0;}
  };

}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
