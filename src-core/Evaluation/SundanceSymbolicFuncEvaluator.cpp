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

#include "SundanceSubtypeEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSet.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;


SymbolicFuncElementEvaluator
::SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr, 
                               const EvalContext& context)
  : SubtypeEvaluator<SymbolicFuncElement>(expr, context),
    mi_(),
    spatialDerivs_(),
    ones_(),
    df_(dynamic_cast<const DiscreteFuncElement*>(expr->evalPt())),
    stringReps_()
{
  
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing symbolic func evaluator for " 
                    << expr->toString());

  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << endl << *(this->sparsity)());

  const ZeroExpr* z 
    = dynamic_cast<const ZeroExpr*>(expr->evalPt());
  
  TEST_FOR_EXCEPTION(z==0 && df_==0, InternalError,
                     "SymbolicFuncElementEvaluator ctor detected an "
                     "evaluation point=" << expr->toString()
                     << " that is neither zero nor a discrete "
                     "function.");

  static Array<string> coordNames;
  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }
  
  int constantCounter = 0;
  int vectorCounter = 0;

  Set<MultiIndex> miSet;
  
  for (int i=0; i<this->sparsity()->numDerivs(); i++) 
    {
      if (this->sparsity()->isSpatialDeriv(i))
        {
          /* evaluate the spatial deriv applied to the evaluation point */
          TEST_FOR_EXCEPTION(z != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "spatial derivative of a zero function. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");

          mi_.append(this->sparsity()->multiIndex(i));
          miSet.put(this->sparsity()->multiIndex(i));
          addVectorIndex(i, vectorCounter);
          spatialDerivs_.append(vectorCounter++);
          int dir = this->sparsity()->multiIndex(i).firstOrderDirection();
          string deriv = "D[" + df_->name() + ", " + coordNames[dir] + "]";
          stringReps_.append(deriv);
        }
      else
        {
          TEST_FOR_EXCEPTION(this->sparsity()->deriv(i).order() > 1,
                             InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "nonzero functional derivative of order greater "
                             "than one. All such derivs should have been "
                             "identified as zero by this point. The bad "
                             "derivative is " << this->sparsity()->deriv(i)
                             << ", and the bad sparsity table is "
                             << *(this->sparsity)());

          if (this->sparsity()->deriv(i).order()==0)
            {
              TEST_FOR_EXCEPTION(z != 0, InternalError,
                             "SymbolicFuncElementEvaluator ctor detected a "
                             "zero-order derivative of a zero function. All "
                             "such expressions should have been "
                             "automatically eliminated by this point.");
              /* value of zeroth functional deriv is a discrete function */
              addVectorIndex(i, vectorCounter);
              spatialDerivs_.append(vectorCounter++);
              mi_.append(MultiIndex());
              miSet.put(MultiIndex());
              stringReps_.append(df_->name());
            }
          else
            {
              /* value of first functional deriv is one */
              addConstantIndex(i, constantCounter);
              ones_.append(constantCounter++);
            }
        }
    }

  if (df_ != 0)
    {
      SUNDANCE_VERB_MEDIUM(tabs << "setting up evaluation for discrete eval pt");
      df_->setupEval(context);
      dfEval_ = dynamic_cast<const DiscreteFuncElementEvaluator*>(df_->evaluator(context).get());
    }
}




void SymbolicFuncElementEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RefCountPtr<EvalVector> >& vectorResults) const 
{
  ///  TimeMonitor timer(symbolicFuncEvalTimer());
  Tabs tabs;
  
  if (verbosity() > VerbSilent)
    {
      cerr << tabs << "SymbolicFuncElementEvaluator::eval: expr=" << expr()->toString() 
           << endl;
      if (verbosity() > VerbLow)
        {
          cerr << tabs << "sparsity = " << endl << *(this->sparsity)() << endl;
        }
    }

  constantResults.resize(ones_.size());
  vectorResults.resize(spatialDerivs_.size());

  /* Evaluate discrete functions if necessary */
  if (df_ != 0 && mi_.size() > 0)
    {
      for (unsigned int i=0; i<mi_.size(); i++)
        {
          vectorResults[i] = mgr.popVector();
          TEST_FOR_EXCEPTION(!vectorResults[i]->isValid(), 
                             InternalError,
                             "invalid evaluation vector allocated in "
                             "SymbolicFuncElementEvaluator::internalEval()");
          vectorResults[i]->setString(stringReps_[i]);
        }
      mgr.evalDiscreteFuncElement(df_, mi_, vectorResults);
    }

  /* Set the known one entries to one */
  for (unsigned int i=0; i<ones_.size(); i++)
    {
      constantResults[ones_[i]] = 1.0;
    }
  

}


