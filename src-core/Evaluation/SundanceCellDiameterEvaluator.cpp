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

#include "SundanceEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

CellDiameterExprEvaluator::CellDiameterExprEvaluator(const CellDiameterExpr* expr, 
                                       const EvalContext& context)
  : SubtypeEvaluator<CellDiameterExpr>(expr, context), 
    stringRep_(expr->toString())
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing cell diameter expr evaluator for " 
                    << expr->toString());
  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *sparsity());

  TEST_FOR_EXCEPTION(sparsity()->numDerivs() > 1, InternalError,
                     "CellDiameterExprEvaluator ctor found a sparsity table "
                     "with more than one entry. The bad sparsity table is "
                     << *sparsity());

  /* 
   * There is only one possible entry in the nozeros table for a
   * cell diameter expression: a zeroth derivative.
   */
  
  for (int i=0; i<sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = sparsity()->deriv(i);

      TEST_FOR_EXCEPTION(d.order()!=0, InternalError,
                         "CellDiameterExprEvaluator ctor found an entry in the "
                         "sparsity superset that is not a zeroth-order derivative. "
                         "The bad entry is " << sparsity()->deriv(i) 
                         << ". The superset is " 
                         << *sparsity());
      addVectorIndex(i, 0);
    }
}



void CellDiameterExprEvaluator::internalEval(const EvalManager& mgr,
                                             Array<double>& constantResults,
                                             Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(cellDiameterEvalTimer());
  Tabs tabs;

  SUNDANCE_VERB_LOW(tabs << "CellDiameterExprEvaluator::eval() expr=" << expr()->toString());

  if (verbosity() > 1)
    {
      cerr << tabs << "sparsity = " << endl << *sparsity() << endl;
    }

  if (sparsity()->numDerivs() > 0)
    {
      vectorResults.resize(1);
      vectorResults[0] = mgr.popVector();
      mgr.evalCellDiameterExpr(expr(), vectorResults[0]);
      if (EvalVector::shadowOps()) vectorResults[0]->setString(stringRep_);
    }

  if (verbosity() > VerbMedium)
    {
      cerr << tabs << "results " << endl;
      sparsity()->print(cerr, vectorResults,
                            constantResults);
    }
}

