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

#ifndef SUNDANCE_EVALMEDIATOR_H
#define SUNDANCE_EVALMEDIATOR_H


#include "SundanceDefs.hpp"
#include "SundanceEvalVector.hpp"
#include "TSFObjectWithVerbosity.hpp"

namespace SundanceCore
{
using namespace SundanceUtils;
class CoordExpr;
class CellDiameterExpr;
class CellVectorExpr;
class MultiIndex; 
class DiscreteFuncElement;
/**
 * Base class for evaluation mediator objects. 
 * Evaluation mediators are responsible
 * for evaluating those expressions whose
 * calculation must be delegated to the framework.
 */
class AbstractEvalMediator 
  : public TSFExtended::ObjectWithVerbosity<AbstractEvalMediator>
{
public:
  /** */
  AbstractEvalMediator(int verb=0);

  /** */
  virtual ~AbstractEvalMediator(){;}

          

  /** Evaluate the given coordinate expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCoordExpr(const CoordExpr* expr,
    RefCountPtr<EvalVector>& vec) const = 0 ;

  /** Evaluate the given discrete function, putting
   * its numerical values in the given EvalVector. */
  virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
    const Array<MultiIndex>& mi,
    Array<RefCountPtr<EvalVector> >& vec) const = 0 ;

  /** Evaluate the given cell diameter expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
    RefCountPtr<EvalVector>& vec) const = 0 ;

  /** Evaluate the given cell vector expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellVectorExpr(const CellVectorExpr* expr,
    RefCountPtr<EvalVector>& vec) const = 0 ;
private:
};
}


#endif
