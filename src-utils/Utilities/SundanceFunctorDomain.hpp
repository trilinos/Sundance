/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_FUNCTORDOMAIN_H
#define SUNDANCE_FUNCTORDOMAIN_H

#include "SundanceDefs.hpp"
#include "Teuchos_StrUtils.hpp"

namespace Sundance
{
  using namespace Teuchos;
  using std::string;

  class FunctorDomain
  {
  public:
    FunctorDomain();

    virtual ~FunctorDomain(){;}

    virtual bool hasLowerBound() const {return false;}

    virtual double lowerBound() const ;

    virtual bool hasUpperBound() const {return false;}

    virtual double upperBound() const ;

    virtual bool hasExcludedPoint() const {return false;}

    virtual double excludedPoint() const ;

    virtual string description() const = 0 ;
  };

  class UnboundedDomain : public FunctorDomain
  {
  public:
    UnboundedDomain();

    string description() const {return "UnboundedDomain()";}
  };


  class PositiveDomain : public FunctorDomain
  {
  public:
    PositiveDomain();

    bool hasLowerBound() const {return true;}
    
    double lowerBound() const {return 0.0;}
    
    string description() const {return "PositiveDomain()";}
  };

  class StrictlyPositiveDomain : public FunctorDomain
  {
  public:
    StrictlyPositiveDomain();

    bool hasLowerBound() const {return true;}
    
    double lowerBound() const {return 0.0;}
    
    bool hasExcludedPoint() const {return true;}
    
    double excludedPoint() const {return 0.0;}

    string description() const {return "StrictlyPositiveDomain()";}

  };


  class BoundedDomain : public FunctorDomain
  {
  public:
    BoundedDomain(const double& lower, const double& upper);

    bool hasLowerBound() const {return true;}
    
    double lowerBound() const {return lower_;}
    
    bool hasUpperBound() const {return true;}
    
    double upperBound() const {return upper_;}

    string description() const {return "BoundedDomain("
	+ Teuchos::toString(lowerBound()) + ", "
	+ Teuchos::toString(upperBound()) + ")";}

  private:
    double lower_;

    double upper_;
  };


  class LowerBoundedDomain : public FunctorDomain
  {
  public:
    LowerBoundedDomain(const double& lower);

     bool hasLowerBound() const {return true;}

     double lowerBound() const {return lower_;}

    string description() const {return "LowerBoundedDomain("
	+ Teuchos::toString(lowerBound()) + ")";}



  private:
    double lower_;
  };

class NonzeroDomain : public FunctorDomain
  {
  public:
    NonzeroDomain();

    bool hasExcludedPoint() const {return true;}
    
    double excludedPoint() const {return 0.0;}

    string description() const {return "NonzeroDomain()";}


  };

  inline std::ostream& operator<<(std::ostream& os, const FunctorDomain& f)
  {
    os << f.description();
    return os;
  }
}



#endif
