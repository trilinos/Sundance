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
/*
 * FunctionalEvaluatorBase.hpp
 *
 *  Created on: Jun 20, 2011
 *      Author: benk
 */

#ifndef FUNCTIONALEVALUATORBASE_HPP_
#define FUNCTIONALEVALUATORBASE_HPP_

#include "SundanceExpr.hpp"

namespace Sundance
{

/** Abstract interface for the FunctionalEvaluator class
 * it contains only empty functions */
class FunctionalEvaluatorBase
{
public:
  /** */
  FunctionalEvaluatorBase() {;}

  /** */
  virtual double evaluate() const { return 0.0; }

  /** */
  virtual Expr evalGradient(double& value) const { return Expr(0.0); }

  /** */
  virtual double fdGradientCheck(double h) const { return 0.0; }
};
}
#endif /* FUNCTIONALEVALUATORBASE_HPP_ */
