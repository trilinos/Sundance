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

#ifndef SUNDANCE_EVENT_DETECTOR_H
#define SUNDANCE_EVENT_DETECTOR_H


#include "SundanceDefs.hpp"

namespace Sundance
{
class Expr;

/** */
class EventDetectorBase
{
public:
  /** */
  EventDetectorBase() {}
  
  /** */
  virtual bool terminateOnDetection() const {return false;}

  /** */
  virtual bool checkForEvent(
    const double& t1, const Expr& u1,
    const double& t2, const Expr& u2) = 0 ;
};

/** */
enum ThresholdEventType {AnyAbove, AllAbove, AnyBelow, AllBelow};

/** */
class ThresholdEventDetector : public EventDetectorBase
{
public:
  /** */
  ThresholdEventDetector(double threshold, ThresholdEventType eventType,
    bool terminateOnDetection=false)
    : threshold_(threshold), eventType_(eventType),
      gotIt_(false), eventTime_(-1.0e300),
      terminateOnDetection_(terminateOnDetection) {}

  /** */
  bool terminateOnDetection() const {return terminateOnDetection_;}

  /** */
  bool checkForEvent(
    const double& t1, const Expr& u1,
    const double& t2, const Expr& u2) ;

  /** */
  double eventTime() const {return eventTime_;}

  /** */
  double foundEvent() const {return gotIt_;}

private:
  double threshold_;
  ThresholdEventType eventType_;
  mutable bool gotIt_;
  mutable double eventTime_;
  bool terminateOnDetection_;
};




}


#endif
