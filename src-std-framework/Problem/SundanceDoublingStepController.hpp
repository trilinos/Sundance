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

#ifndef SUNDANCE_DOUBLING_STEP_CONTROLLER_H
#define SUNDANCE_DOUBLING_STEP_CONTROLLER_H

#include "SundanceFieldWriter.hpp"
#include "SundanceFieldWriterFactory.hpp"
#include "SundanceTransientStepProblem.hpp"
#include "SundanceExprComparison.hpp"
#include "SundanceEventDetector.hpp"



namespace Sundance
{

/** */
class StepControlParameters
{
public:
  /** */
  StepControlParameters()
    {initDefaults();}

  double tStart_;

  double tStop_;

  double tau_;

  double initialStepsize_;

  double minStepsizeFactor_;

  double maxStepsizeFactor_;

  double stepsizeReductionSafetyFactor_;

  int maxSteps_;

  int verbosity_;

  int stepOrder_;

private:
  void initDefaults()
    {
      tStart_ = 0.0;
      tStop_ = 0.0;
      tau_ = 1.0e-6;
      initialStepsize_ = 0.01;
      minStepsizeFactor_ = 0.01;
      maxStepsizeFactor_ = 10.0;
      stepsizeReductionSafetyFactor_ = 0.9;
      maxSteps_ = 100000;
      verbosity_ = 0;
      stepOrder_ = 2;
    }
};

/** */
class StepHookBase
{
public:
  virtual void call(const double& tCur, const Expr& uCur) const = 0 ;
};

/** */
class OutputControlParameters
{
public:
  /** */
  OutputControlParameters(
    const FieldWriterFactory& wf,
    const string& filename,
    const double& writeInterval,
    int verb=0)
    : 
    writeInterval_(writeInterval),
    wf_(wf),
    filename_(filename),
    verbosity_(verb)
    {}
  
  double writeInterval_;
  FieldWriterFactory wf_;
  string filename_;
  int verbosity_;
};


/** */
class DoublingStepController
{
public:
  /** */
  DoublingStepController(
    const TransientStepProblem& prob,
    const NonlinearSolver<double>& solver,
    const StepControlParameters& stepControl,
    const OutputControlParameters& outputControl,
    const RCP<ExprComparisonBase>& compare)
    : prob_(prob), 
      stepControl_(stepControl), 
      outputControl_(outputControl),
      solver_(solver),
      compare_(compare),
      eventHandler_()
    {}

  /** */
  void setEventHandler(RCP<EventDetectorBase> e)
    {eventHandler_ = e;}

  /** */
  void setStepHook(RCP<StepHookBase> h)
    {stepHook_ = h;}
      

  /** */
  bool run() const ;

  /** */
  void write(int index, double t, const Expr& u) const ;

private:
  TransientStepProblem prob_;
  StepControlParameters stepControl_;
  OutputControlParameters outputControl_;
  NonlinearSolver<double> solver_;
  RCP<ExprComparisonBase> compare_;
  RCP<EventDetectorBase> eventHandler_;
  RCP<StepHookBase> stepHook_;
};


}


#endif
