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

#include "SundanceTransientStepProblem.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

namespace Sundance
{

TransientStepProblem::TransientStepProblem(
    const NonlinearProblem& stepProb, 
    const Expr& tPrev, const Expr& uPrev, 
    const Expr& tNext, const Expr& uNext,
    const Expr& dt,
    int verbosity
    )
  : prob_(stepProb),
    tPrev_(tPrev),
    tNext_(tNext),
    uPrev_(uPrev),
    uNext_(uNext),
    dt_(dt),
    verb_(verbosity)
{
}

bool TransientStepProblem::step(
  double tCur, const Expr& uCur,
  double tNext, Expr uNext,
  const NonlinearSolver<double>& solver) const
{
  Tabs tab;
  
  PLAYA_MSG1(verb_, tab << "step from t=" << tCur << " to " << tNext);
  tPrev_.setParameterValue(tCur);
  tNext_.setParameterValue(tNext);
  dt_.setParameterValue(tNext-tCur);

  PLAYA_MSG2(verb_, tab << "updating uPrev");
  updateDiscreteFunction(uCur, uPrev_);

  PLAYA_MSG2(verb_, tab << "updating uNext");
  updateDiscreteFunction(uNext, uNext_);

  PLAYA_MSG2(verb_, tab << "Solving NLP");
  SolverState<double> state = prob_.solve(solver);

  PLAYA_MSG2(verb_, tab << "Updating solution");
  updateDiscreteFunction(uNext_, uNext);

  return state.finalState() == SolveConverged;
}

}
