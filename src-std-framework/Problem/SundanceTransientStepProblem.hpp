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

#ifndef SUNDANCE_TRANSIENT_STEP_PROBLEM_H
#define SUNDANCE_TRANSIENT_STEP_PROBLEM_H


#include "SundanceDefs.hpp"
#include "SundanceNonlinearProblem.hpp"

namespace Sundance
{
using Playa::NonlinearSolver;

/** 
 *
 */
class TransientStepProblem
{
public:
  /** */
  TransientStepProblem(
    const NonlinearProblem& stepProb, 
    const Expr& tPrev, const Expr& uPrev, 
    const Expr& tNext, const Expr& uNext,
    const Expr& dt,
    int verbosity=0
    );

  /** */
  bool step(double tCur, const Expr& uCur,
    double tNext, Expr uNext,
    const NonlinearSolver<double>& solver) const ;


  /** */
  Expr uCur() const {return uPrev_;}
  /** */
  Expr uSoln() const {return uNext_;}


private:
  NonlinearProblem prob_;
  mutable Expr tPrev_;
  mutable Expr tNext_;
  mutable Expr uPrev_;
  mutable Expr uNext_;
  mutable Expr dt_;
  int verb_;
};

}


#endif
